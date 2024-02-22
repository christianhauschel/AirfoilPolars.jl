"""
    AirfoilPolars

This module provides tools to handle airfoil polars.

The module provides the following functions:

- `Polar`: constructor for the `Polar` type
- `extrapolate`: extrapolate the polars
- `blend`: blend two polars
- `smooth`: smooth the polars
- `copy`: copy a polar
- `save`: save a polar to a file
- `plot`: plot a polar
- `generate_name`: generate a name for the polar
- `correction3D`: correct the 2D polars to 3D
- `correction_Mach`: correct the polars for Mach number
"""
module AirfoilPolars

include("polar.jl")
export Polar, extrapolate, blend, smooth, copy, save, plot, generate_name
export correction3D, correction_Mach

include("io.jl")
export get_filename_ext, load

include("solve.jl")
export init, solve

end
