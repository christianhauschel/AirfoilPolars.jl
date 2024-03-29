using LinearAlgebra
using Interpolations
using Statistics
using FLOWMath
import Dierckx
using CSV
using PyFormattedStrings
using DataFrames
using LsqFit
using Plots, LaTeXStrings

"""
    Polar

A struct to store airfoil polars.

# Fields
- `Re::Float64`: Reynolds number
- `alpha::Vector{Float64}`: angle of attack [deg]
- `cl::Vector{Float64}`: lift coefficient
- `cd::Vector{Float64}`: drag coefficient
- `cm::Vector{Float64}`: moment coefficient
- `M::Float64`: Mach number
- `n_crit::Float64`: critical amplification factor
- `xtrip::Tuple{Float64,Float64}`: transition location
- `name_airfoil::String`: name of the airfoil
"""
mutable struct Polar
    Re::Float64
    alpha::Vector{Float64}
    cl::Vector{Float64}
    cd::Vector{Float64}
    cm::Vector{Float64}
    M::Float64
    n_crit::Float64
    xtrip::Tuple{Float64,Float64}
    name_airfoil::String
end

"""
    Polar(Re, alpha, cl, cd, cm)

Create a new Polar object with no Xfoil properties.

# Arguments
- `Re::Float64`: Reynolds number
- `alpha::Vector{Float64}`: angle of attack [deg]
- `cl::Vector{Float64}`: lift coefficient
- `cd::Vector{Float64}`: drag coefficient
- `cm::Vector{Float64}`: moment coefficient
"""
function Polar(Re, alpha, cl, cd, cm)
    return Polar(Re, alpha, cl, cd, cm, NaN64, NaN64, (NaN64, NaN64), "NONAME")
end

"""
    smooth1D(x, y; order=3, smoothing=0.01)

Smooths a 1D curve using splines.
"""
function _smooth1D(x, y; order=3, smoothing=0.01)
    crv = Dierckx.Spline1D(x, y; k=order, s=smoothing)
    return crv(x)
end

function Base.copy(p::Polar)
    return Polar(p.Re, copy(p.alpha), copy(p.cl), copy(p.cd), copy(p.cm), copy(p.M), copy(p.n_crit), copy(p.xtrip), p.name_airfoil)
end

"""
    smooth(
        p::Polar; 
        order_cl=3, order_cd=3, order_cm=3, 
        smoothing_cl=0.01, smoothing_cd=0.0001, smoothing_cm=0.005
    )

Smooths a polar using splines.
"""
function smooth(
    p::Polar;
    order_cl=3,
    order_cd=3,
    order_cm=3,
    smoothing_cl=0.1,
    smoothing_cd=0.01,
    smoothing_cm=0.005
)
    cl = _smooth1D(p.alpha, p.cl; order=order_cl, smoothing=smoothing_cl)
    cd = _smooth1D(p.alpha, p.cd; order=order_cd, smoothing=smoothing_cd)
    cm = _smooth1D(p.alpha, p.cm; order=order_cm, smoothing=smoothing_cm)

    return Polar(copy(p.Re), copy(p.alpha), cl, cd, cm, copy(p.M), copy(p.n_crit), (copy(p.xtrip[1]), copy(p.xtrip[2])), p.name_airfoil)
end

"""
    blend(p1::Polar, p2::Polar, weight2::Float64)

Blends two polars using a weight.

# Arguments
- `p1::Polar`: first polar
- `p2::Polar`: second polar
- `weight2::Float64`: weight for the second polar
"""
function blend(p1::Polar, p2::Polar, weight2::Float64)
    alpha = union(p1.alpha, p2.alpha)

    min_alpha = maximum([minimum(p1.alpha), minimum(p2.alpha)])
    max_alpha = minimum([maximum(p1.alpha), maximum(p2.alpha)])
    alpha = alpha[(alpha.>=min_alpha).&(alpha.<=max_alpha)]

    cl1 = akima(p1.alpha, p1.cl, alpha)
    cl2 = akima(p2.alpha, p2.cl, alpha)
    cd1 = akima(p1.alpha, p1.cd, alpha)
    cd2 = akima(p2.alpha, p2.cd, alpha)
    cm1 = akima(p1.alpha, p1.cm, alpha)
    cm2 = akima(p2.alpha, p2.cm, alpha)

    Re = p1.Re + weight2 * (p2.Re - p1.Re)
    cl = cl1 + weight2 * (cl2 - cl1)
    cd = cd1 + weight2 * (cd2 - cd1)
    cm = cm1 + weight2 * (cm2 - cm1)

    return Polar(Re, alpha, cl, cd, cm, mean([p1.M, p2.M]), mean([p1.n_crit, p2.n_crit]), mean.([p1.xtrip, p2.xtrip]))
end

"""
    _correction3D(
        alpha::Vector{Float64}, cl_2d::Vector{Float64}, cd_2d::Vector{Float64}, 
        r::Float64, c::Float64, u_inf::Float64, Ω::Float64; 
        alpha_max_corr::Float64=30.0, alpha_linear_min::Float64=-5.0, alpha_linear_max::Float64=5.0
    )

Applies 3-D corrections for rotating sections from the 2-D data using the Du and Selig model.

# Parameters
- `alpha::Vector{Float64}`: angle of attack [deg]
- `cl_2d::Vector{Float64}`: 2-D lift coefficient
- `cd_2d::Vector{Float64}`: 2-D drag coefficient
- `r::Float64`: local radial location
- `c::Float64`: local chord 
- `u_inf::Float64`: local sectional airfoil freestream velocity ≈ √(Vwind^2 + (Ωr)^2)
- `Ω::Float64`: rotation speed [rad/s]
- `alpha_max_corr::Float64`: maximum angle of attack to apply full correction
- `alpha_linear_min::Float64`: angle of attack where linear portion of lift curve slope begins
- `alpha_linear_max::Float64`: angle of attack where linear portion of lift curve slope ends
"""
function _correction3D(
    alpha::Vector{Float64}, cl_2d::Vector{Float64}, cd_2d::Vector{Float64},
    r::Float64, c::Float64, u_inf::Float64, Ω::Float64;
    alpha_max_corr::Float64=30.0, alpha_linear_min::Float64=-0.0, alpha_linear_max::Float64=5.0
)
    # Convert units for convenience
    alpha = deg2rad.(alpha)
    alpha_max_corr = deg2rad(alpha_max_corr)
    alpha_linear_min = deg2rad(alpha_linear_min)
    alpha_linear_max = deg2rad(alpha_linear_max)

    # Parameters in Du-Selig model
    a = 1
    b = 1
    d = 1

    λ = Ω * r / u_inf
    expon = d / λ #d / Λ / r_over_R
    expon_d = expon / 2

    n_alpha = length(alpha)

    # Find linear region
    idx = (alpha .>= alpha_linear_min) .& (alpha .<= alpha_linear_max)

    # fit line to region using Least Squares
    model(t, p) = p[1] * t .+ p[2]

    fit = curve_fit(model, alpha[idx], cl_2d[idx], [2π, 0.0])
    m = fit.param[1]
    alpha0 = -fit.param[2] / m

    chord_over_r = c / r

    # Correction factor
    fcl = 1.0 / m * (1.6 * chord_over_r / 0.1267 * (a - chord_over_r .^ expon) ./ (b + chord_over_r .^ expon) - 1)
    fcd = 1.0 / m * (1.6 * chord_over_r / 0.1267 * (a - chord_over_r .^ expon_d) ./ (b + chord_over_r .^ expon_d) - 1)

    # Adjustment
    amax = atan(1 / 0.12) - deg2rad(5)
    adj = zeros(n_alpha)
    for i = 1:n_alpha
        if abs(alpha[i]) >= amax
            adj[i] = 0.0
        elseif abs(alpha[i]) > alpha_max_corr
            adj[i] = ((amax - abs(alpha[i])) / (amax - alpha_max_corr))^2
        else
            adj[i] = 1.0
        end
    end

    # Du-Selig correction for lift
    cl_linear = m * (alpha .- alpha0)
    Δcl = fcl .* (cl_linear - cl_2d) .* adj
    cl_3d = cl_2d .+ Δcl

    # Du-Selig correction for drag
    cd0 = linear(alpha, cd_2d, 0.0)
    Δcd = fcd .* (cd0 .- cd_2d)
    cd_3d = cd_2d .- Δcd

    return cl_3d, cd_3d
end


"""
    correction3D(
        polar::Polar, 
        r::Float64, c::Float64, u_inf::Float64, Ω::Float64; 
        alpha_max_corr::Float64=30.0, alpha_linear_min::Float64=-5.0, alpha_linear_max::Float64=5.0
    )

Applies 3-D corrections for rotating sections from the 2-D data using the Du and Selig model.


# Parameters 
- `polar::Polar`: 2-D polar
- `r::Float64`: local radius
- `c::Float64`: local chord 
- `u_inf::Float64`: local airfoil freestream velocity ≈ √(Vwind^2 + (Ωr)^2)
- `Ω::Float64`: rotation speed [rad/s]
- `alpha_max_corr::Float64`: maximum angle of attack to apply full correction
- `alpha_linear_min::Float64`: angle of attack where linear portion of lift curve slope begins
- `alpha_linear_max::Float64`: angle of attack where linear portion of lift curve slope ends

# Returns
- `polar::Polar`: A new Polar object corrected for 3-D effects

# References
1. Z. Du and M. Selig, ‘A 3-D stall-delay model for horizontal axis wind turbine performance prediction’, in 1998 ASME Wind Energy Symposium, Reno, NV, U.S.A.: American Institute of Aeronautics and Astronautics, Jan. 1998. doi: 10.2514/6.1998-21.
"""
function correction3D(
    polar::Polar, r::Float64, c::Float64, u_inf::Float64, Ω::Float64;
    alpha_max_corr::Float64=30.0, alpha_linear_min::Float64=-5.0, alpha_linear_max::Float64=5.0
)

    cl_3d, cd_3d = _correction3D(
        polar.alpha, polar.cl, polar.cd,
        r, c, u_inf, Ω;
        alpha_max_corr=alpha_max_corr, alpha_linear_min=alpha_linear_min, alpha_linear_max=alpha_linear_max
    )

    return Polar(polar.Re, polar.alpha, cl_3d, cd_3d, polar.cm, polar.M, polar.n_crit, polar.xtrip, polar.name_airfoil * "_3D")
end

function _correction_Mach(cl, Ma_new)
    cl = cl ./ sqrt(1 - Ma_new^2)
    return cl
end

"""
    correction_Mach(polar::Polar, Ma_new)

Corrects a polar for a new Mach number using Prandtl-Glauert compressibility correction:

``c_l = \\frac{c_{l, M=0}}{\\sqrt{1 - M^2}}``
"""
function correction_Mach(polar::Polar, Ma_new)
    cl = _correction_Mach(polar.cl, Ma_new)
    return Polar(polar.Re, polar.alpha, cl, polar.cd, polar.cm, Ma_new, polar.n_crit, polar.xtrip, polar.name_airfoil)
end

function _Viterna(alpha::Vector, cl_adj, cd_max, A, B)
    alpha = max.(alpha, 0.0001)
    cl = cd_max / 2 * sin.(2 .* alpha) + A * cos.(alpha) .^ 2 ./ sin.(alpha)
    cl *= cl_adj
    cd = cd_max * sin.(alpha) .^ 2 + B * cos.(alpha)
    return cl, cd
end

function _CMCoeff(polar::Polar, cl_high::Float64, cd_high::Float64, cm_high::Float64)
    found_zero_lift = false
    cm0 = 0.0

    for i in 1:(length(polar.cm)-1)
        if abs(polar.alpha[i]) < 20.0 && polar.cl[i] <= 0 && polar.cl[i+1] >= 0
            p = -polar.cl[i] / (polar.cl[i+1] - polar.cl[i])
            cm0 = polar.cm[i] + p * (polar.cm[i+1] - polar.cm[i])
            found_zero_lift = true
            break
        end
    end

    if !found_zero_lift
        p = -polar.cl[1] / (polar.cl[2] - polar.cl[1])
        cm0 = polar.cm[1] + p * (polar.cm[2] - polar.cm[1])
    end

    alpha_high = deg2rad(polar.alpha[end])
    XM = (-cm_high + cm0) / (cl_high * cos(alpha_high) + cd_high * sin(alpha_high))
    cmCoef = (XM - 0.25) / tan(alpha_high - pi / 2)

    return cmCoef, cm0
end


function _getCM(p::Polar, i::Int, cm0, cmCoef, alpha::Vector, cl_ext::Vector, cd_ext::Vector, alpha_low_deg, alpha_high_deg)
    cm_new = 0.0

    if alpha[i] >= alpha_low_deg && alpha[i] <= alpha_high_deg
        return
    end

    if alpha[i] > -165 && alpha[i] < 165
        if abs(alpha[i]) < 0.01
            cm_new = cm0
        else
            if alpha[i] > 0
                x = cmCoef * tan(deg2rad(alpha[i]) - pi / 2) + 0.25
                cm_new = cm0 - x * (cl_ext[i] * cos(deg2rad(alpha[i])) + cd_ext[i] * sin(deg2rad(alpha[i])))
            else
                x = cmCoef * tan(-deg2rad(alpha[i]) - pi / 2) + 0.25
                cm_new = -(cm0 - x * (-cl_ext[i] * cos(-deg2rad(alpha[i])) + cd_ext[i] * sin(-deg2rad(alpha[i]))))
            end
        end
    else
        if alpha[i] == 165
            cm_new = -0.4
        elseif alpha[i] == 170
            cm_new = -0.5
        elseif alpha[i] == 175
            cm_new = -0.25
        elseif alpha[i] == 180
            cm_new = 0.0
        elseif alpha[i] == -165
            cm_new = 0.35
        elseif alpha[i] == -170
            cm_new = 0.4
        elseif alpha[i] == -175
            cm_new = 0.2
        elseif alpha[i] == -180
            cm_new = 0.0
        else
            println("Angle encountered for which there is no CM table value (near +/-180 deg). Program will stop.")
        end
    end

    return cm_new
end


"""
    extrapolate(
        p::Polar; 
        cd_max::Union{Nothing,Float64}=nothing, 
        AR::Union{Nothing,Float64}=nothing, 
        cd_min=0.001, 
        n_alpha::Int=15
    )

Extrapolate a polar to +/- 180 degrees. This function is based on the Viterna method.

# Arguments 
- `p::Polar`: polar to extrapolate
- `cd_max::Union{Nothing,Float64}`: maximum drag coefficient
- `AR::Union{Nothing,Float64}`: aspect ratio
- `cd_min::Float64`: minimum drag coefficient
- `n_alpha::Int`: number of points to use for extrapolation
"""
function extrapolate(
    p::Polar;
    cd_max::Union{Nothing,Float64}=nothing,
    AR::Union{Nothing,Float64}=nothing,
    cd_min=0.001,
    n_alpha::Int=15
)
    if cd_min < 0
        error("cdmin cannot be < 0")
    end

    # Lift coefficient adjustment to account for asymmetry
    cl_adj = 0.7

    # Estimate CD max
    if cd_max === nothing
        if AR !== nothing
            cd_max = 1.11 + 0.018 * AR
        else
            error("Please provide AR or cdmax!")
        end
    end
    cd_max = max(maximum(p.cd), cd_max)

    # Extract matching info from ends
    alpha_high = deg2rad(p.alpha[end])
    cl_high = p.cl[end]
    cd_high = p.cd[end]
    cm_high = p.cm[end]

    alpha_low = deg2rad(p.alpha[1])
    cl_low = p.cl[1]
    cd_low = p.cd[1]

    if alpha_high > pi / 2 || alpha_low < -pi / 2
        error("alpha[-1] > pi/2 or alpha[0] < -pi/2")
    end

    # Parameters used in model
    sa = sin(alpha_high)
    ca = cos(alpha_high)
    A = (cl_high - cd_max * sa * ca) * sa / ca^2
    B = (cd_high - cd_max * sa^2) / ca

    # alpha_high <-> 90
    alpha1 = Vector(range(alpha_high, stop=pi / 2, length=n_alpha))[2:end]
    cl1, cd1 = _Viterna(alpha1, 1.0, cd_max, A, B)

    # 90 <-> 180-alpha_high
    alpha2 = Vector(range(pi / 2, stop=pi - alpha_high, length=n_alpha))[2:end]
    cl2, cd2 = _Viterna(pi .- alpha2, -cl_adj, cd_max, A, B)

    # 180-alpha_high <-> 180
    alpha3 = Vector(range(pi - alpha_high, stop=pi, length=n_alpha))[2:end]
    cl3, cd3 = _Viterna(pi .- alpha3, 1.0, cd_max, A, B)
    cl3 .= (alpha3 .- pi) ./ alpha_high .* cl_high .* cl_adj  # override with linear variation

    if alpha_low <= -alpha_high
        alpha4 = Float64[]
        cl4 = Float64[]
        cd4 = Float64[]
        alpha5max = alpha_low
    else
        # -alpha_high <-> alpha_low
        alpha4 = Vector(range(-alpha_high, stop=alpha_low, length=n_alpha))[2:end-1]  # also remove last element for concatenation for this case
        cl4 = -cl_high .* cl_adj .+ (alpha4 .+ alpha_high) ./ (alpha_low .+ alpha_high) .* (cl_low .+ cl_high .* cl_adj)
        cd4 = cd_low .+ (alpha4 .- alpha_low) ./ (-alpha_high .- alpha_low) .* (cd_high .- cd_low)
        alpha5max = -alpha_high
    end

    # -90 <-> -alpha_high
    alpha5 = Vector(range(-pi / 2, stop=alpha5max, length=n_alpha))[2:end]
    cl5, cd5 = _Viterna(-alpha5, -cl_adj, cd_max, A, B)

    # -180+alpha_high <-> -90
    alpha6 = Vector(range(-pi + alpha_high, stop=-pi / 2, length=n_alpha))[2:end]
    cl6, cd6 = _Viterna(alpha6 .+ pi, cl_adj, cd_max, A, B)

    # -180 <-> -180 + alpha_high
    alpha7 = Vector(range(-pi, stop=-pi + alpha_high, length=n_alpha))
    cl7, cd7 = _Viterna(alpha7 .+ pi, 1.0, cd_max, A, B)
    cl7 .= (alpha7 .+ pi) ./ alpha_high .* cl_high .* cl_adj  # linear variation


    # Concatenating alpha, cl, cd arrays
    alpha = vcat(alpha7, alpha6, alpha5, alpha4, deg2rad.(p.alpha), alpha1, alpha2, alpha3)
    cl = vcat(cl7, cl6, cl5, cl4, p.cl, cl1, cl2, cl3)
    cd = vcat(cd7, cd6, cd5, cd4, p.cd, cd1, cd2, cd3)
    cd = max.(cd, cd_min)  # Don't allow negative drag coefficients

    # ---------------------------------
    # CM extrapolation
    # ---------------------------------

    # Setup alpha and cm to be used in extrapolation
    eps = 1e-13
    cm1_alpha = floor(p.alpha[1] / 10.0 - eps) * 10.0
    cm2_alpha = ceil(p.alpha[end] / 10.0 + eps) * 10.0
    alpha_num = abs(Int((-180.0 - cm1_alpha) / 10.0 - 1))
    alpha_cm1 = Vector(range(-180.0, stop=cm1_alpha, length=alpha_num))
    alpha_cm2 = Vector(range(cm2_alpha, stop=180.0, length=Int((180.0 - cm2_alpha) / 10.0 + 1)))
    alpha_cm = vcat(alpha_cm1, p.alpha, alpha_cm2)  # Specific alpha values are needed for cm function to work
    cm1 = zeros(length(alpha_cm1))
    cm2 = zeros(length(alpha_cm2))
    cm_ext = vcat(cm1, p.cm, cm2)

    if count(!iszero, p.cm) > 0
        cmCoef, cm0 = _CMCoeff(p, cl_high, cd_high, cm_high)  # Get cm coefficient
        cl_cm = akima(rad2deg.(alpha), cl, alpha_cm)  # Get cl for applicable alphas
        cd_cm = akima(rad2deg.(alpha), cd, alpha_cm)  # Get cd for applicable alphas
        alpha_low_deg = p.alpha[1]
        alpha_high_deg = p.alpha[end]

        for i in 1:length(alpha_cm)

            cm_new = _getCM(p, i, cm0, cmCoef, alpha_cm, cl_cm, cd_cm, alpha_low_deg, alpha_high_deg)
            if cm_new !== nothing
                cm_ext[i] = cm_new
            end
        end
    end
    cm = akima(alpha_cm, cm_ext, rad2deg.(alpha))

    return Polar(p.Re, rad2deg.(alpha), cl, cd, cm, p.M, p.n_crit, p.xtrip, p.name_airfoil)
end

"""
    generate_name(p::Polar; fname_extra="")

Generates a name for a polar based on its properties.

# Arguments
- `p::Polar`: polar to generate name for
- `fname_extra::String`: extra string to append to the name
"""
function generate_name(p::Polar; fname_extra="")
    if (p.alpha[1] == -180 && p.alpha[end] == 180) | (p.alpha[1] == -π && p.alpha[end] == π)
        str_extra = "_360"
    else
        str_extra = ""
    end
    if p.name_airfoil == ""
        name_af = "NONAME"
    else
        name_af = p.name_airfoil
    end
    if fname_extra != ""
        fname_extra = "_" * fname_extra
    end

    return f"{name_af}_Re{p.Re/1e6:0.3f}_M{p.M:0.2f}_N{p.n_crit:0.1f}{str_extra}{fname_extra}"
end

"""
    plot(polars::Vector{Polar}; fname=nothing, dpi=300, legend=true)

Plots a vector of polars.
"""
function plot(polars::Vector{Polar}; fname=nothing, dpi=300, legend=true)
    plot_layout = Plots.plot(layout=(3, 1), size=(600, 800), dpi=dpi)

    names_polar = generate_name.(polars)
    names_polar = reshape(names_polar, 1, length(names_polar))

    # re = [f"{p.Re:0.2e}" for p in polars]
    # re = reshape(re, 1, length(re))
    cls = [p.cl for p in polars]
    cds = [p.cd for p in polars]
    cms = [p.cm for p in polars]
    alphas = [p.alpha for p in polars]

    plot!(plot_layout[1], alphas, cls, label=names_polar, ylabel=L"$c_l$", legend=legend)
    plot!(plot_layout[2], alphas, cds, ylabel=L"$c_d$", legend=false)
    plot!(plot_layout[3], alphas, cms, ylabel=L"$c_m$", legend=false, xlabel=L"$\alpha$ [deg]")

    # plot!(plot_layout[1], legend=:bottomright, legendtitle="Reynolds")

    display(plot_layout)

    if fname !== nothing
        savefig(plot_layout, fname)
    end
end