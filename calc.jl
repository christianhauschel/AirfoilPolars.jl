using AirfoilPolars
import AirfoilFast as AF
using PrettyPrint
using PyFormattedStrings

# ==============================================================================
# Settings
# ==============================================================================

dir_data = "/mnt/a/Code/10_drones/Propeller/data"

dir_out = join(dir_data, "polars")
af = AF.Airfoil("data/airfoils/NACA0012.csv")

cd_max = 1.5
n_re = 2
re = Vector(LinRange(20e3, 100e3, n_re))
alpha = Vector(0.0:0.5:4)
n_re = length(re)
mach = 0.1

# ==============================================================================

polars = Vector{Polar}(undef, n_re)

init(af)

section("Solve")
for i in 1:n_re
    subsubsection(f"Re: {re[i]/1000:7.1f} k")
    polars[i] = solve(alpha, re[i]; mach=mach)
end

polars_smooth = smooth.(polars)



polars_ext = extrapolate.(polars_smooth; cd_max = cd_max )

names_polar = generate_name.(polars_ext)

plot([polars; polars_ext]; fname="plot.png")

save.(polars_ext, [joinpath(dir_out, name*".csv") for name in names_polar])