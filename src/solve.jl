import Xfoil
using AirfoilFast
using .AirfoilPolars: Polar
using FLOWMath: linear

name_af = ""

function init(af::Airfoil)
    Xfoil.set_coordinates(af.x, af.y)
    global name_af
    name_af = af.name
end

function solve(
    alpha::Vector,
    re::Number;
    n_iter::Int = 100,
    n_crit::Number = 9.0,
    xtrip::Tuple = (1.0, 1.0),
    mach::Number = 0.0,
    make_nonconverged_nan::Bool = false,    
    interpolate_nonconverged::Bool = false,
)
    n_alpha = length(alpha)
    cl = zeros(n_alpha)
    cd = zeros(n_alpha)
    cdp = zeros(n_alpha)
    cm = zeros(n_alpha)
    converged = zeros(n_alpha)
    for i = 1:n_alpha
        cl[i], cd[i], cdp[i], cm[i], converged[i] = Xfoil.solve_alpha(
            alpha[i],
            re;
            iter = n_iter,
            reinit = i == 1,
            ncrit = n_crit,
            xtrip = xtrip,
            mach = mach,
        )
    end
    global name_af

    if make_nonconverged_nan
        cl[converged .== 0] .= NaN
        cd[converged .== 0] .= NaN
        cm[converged .== 0] .= NaN
    end

    if interpolate_nonconverged
        cl[converged .== 0] .= NaN
        cd[converged .== 0] .= NaN
        cm[converged .== 0] .= NaN
        cl = interpolate_nan(alpha, cl)
        cd = interpolate_nan(alpha, cd)
        cm = interpolate_nan(alpha, cm)
    end

    return Polar(re, alpha, cl, cd, cm, mach, n_crit, xtrip, name_af)
end

function interpolate_nan(x::Vector, y::Vector)
    # select only values that are not NaN 
    y_valid = y[.!isnan.(y)]
    x_valid = x[.!isnan.(y)]

    # interpolate  
    y_new = akima(x_valid, y_valid, x)

    return y_new
end