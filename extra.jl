using AirfoilPolars
using DataFrames, CSV

df = CSV.read("data/polars/NACA4408_Re0.020_M0.00_N9.0.csv", DataFrame)

p = Polar(df.re[1], df.alpha, df.cl, df.cd, df.cm)
p_smooth = smooth(p)
p_ext = extrapolate(p_smooth; cd_max=1.5)

plot([p_ext, p])