using AirfoilPolars
import AirfoilFast as AF
using PrettyPrint
using PyFormattedStrings

# ==============================================================================
# Settings
# ==============================================================================

dir_data = "/home/hacs/projects/neaptide/Propeller/data"

dir_out = joinpath(dir_data, "polars")
af = AF.Airfoil(joinpath(dir_data, "airfoils/NACA4412.dat"))

cd_max = 1.5
# n_re = 2
# re = Vector(LinRange(20e3, 100e3, n_re))
# re = [100e3]
re = [1e6]
n_re = length(re)
alpha = Vector(-14:0.25:18.75)
n_re = length(re)
mach = 0.0
n_crit = 9
n_iter = 100

# ==============================================================================

polars = Vector{Polar}(undef, n_re)


init(af)

section("Solve")
for i in 1:n_re
    subsubsection(f"Re: {re[i]/1000:7.1f} k")
    polars[i] = solve(alpha, re[i]; mach=mach, n_crit = n_crit, make_nonconverged_nan=false, interpolate_nonconverged=false, n_iter=n_iter)
end

plot(polars)

polars_smooth = smooth.(polars; smoothing_cl=0.1)
plot([polars; polars_smooth]; legend=false)

polars_ext = extrapolate.(polars_smooth; cd_max = cd_max )

Ω = 6000 / 60 * 2π
Vwind = 10
c = 0.013
r = 0.127 * 0.75
tsr = Ω * r / Vwind
u_inf = sqrt((Ω * r)^2+Vwind^2)

@time polars_3d = correction3D.(polars_smooth, r, c, u_inf, Ω);

plot([polars_smooth[1], polars_3d[1]])

polars_Ma = correction_Mach.(polars_smooth, 0.4)

plot([polars_smooth[1], polars_Ma[1]]; fname="docs/img/mach.svg")


plot([polars_smooth[1]; polars_3d[1]])

names_polar = generate_name.(polars)
names_polar_ext = generate_name.(polars_ext)

# plot([polars; polars_ext]; fname="plot.png")

polars_ext_3d = correction3D.(polars_ext, r, c, u_inf, Ω);
# polars_ext_3d = correction3D.(polars_ext, r, c, u_inf, Ω);


plot([polars_ext[1]; polars_ext_3d[1]]; legend=false, fname="docs/img/plot.svg")

# save.(polars, [joinpath(dir_out, name*".csv") for name in names_polar])
# save.(polars_ext, [joinpath(dir_out, name*".csv") for name in names_polar_ext])

