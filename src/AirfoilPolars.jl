module AirfoilPolars

include("polar.jl")
export Polar, extrapolate, blend, smooth, copy, save, plot, generate_name

include("solve.jl")
export init, solve

end
