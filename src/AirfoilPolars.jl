module AirfoilPolars

include("polar.jl")
export Polar, extrapolate, blend, smooth, copy, save, plot, generate_name
export correction3D, correction_Mach

include("io.jl")
export get_filename_ext, load

include("solve.jl")
export init, solve

end
