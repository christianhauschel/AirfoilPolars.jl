using AirfoilPolars
import AirfoilFast as AF
using PrettyPrint
using PyFormattedStrings

# ==============================================================================
# Settings
# ==============================================================================

dir_data = "/mnt/a/Code/10_drones/Propeller/data"

dir_out = joinpath(dir_data, "polars")
af = AF.Airfoil(joinpath(dir_data, "airfoils/NACA4412.dat"))

cd_max = 1.5
# n_re = 2
# re = Vector(LinRange(20e3, 100e3, n_re))
# re = [100e3]
re = [50e3]
n_re = length(re)
alpha = Vector(-10:0.25:15)
n_re = length(re)
mach = 0.0
n_crit = 5
n_iter = 100

# ==============================================================================

polars = Vector{Polar}(undef, n_re)

init(af)

section("Solve")
for i in 1:n_re
    subsubsection(f"Re: {re[i]/1000:7.1f} k")
    polars[i] = solve(alpha, re[i]; mach=mach, n_crit = n_crit, make_nonconverged_nan=false, interpolate_nonconverged=false, n_iter=n_iter)
end

polars_smooth = smooth.(polars; smoothing_cm=0.00025)

plot([polars; polars_smooth])

polars_ext = extrapolate.(polars_smooth; cd_max = cd_max )

@time polars_3d = correction3D.(polars_smooth, 0.22, 0.53, 4.24);

fig = plot([polars_smooth[1]; polars_3d[1]])

names_polar = generate_name.(polars)
names_polar_ext = generate_name.(polars_ext)

plot([polars; polars_ext]; fname="plot.png")

polars_ext_3d = correction3D.(polars_ext, 0.22, 0.53, 4.24);
polars_ext_3d = correction3D.(polars_ext, 0.12037037037037036, 1.1626923076923077, 8.4823);


plot([polars_ext[1]; polars_ext_3d[1]])

# save.(polars, [joinpath(dir_out, name*".csv") for name in names_polar])
# save.(polars_ext, [joinpath(dir_out, name*".csv") for name in names_polar_ext])

