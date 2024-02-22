using Revise, Documenter, AirfoilPolars

makedocs(
    sitename="AirfoilPolars.jl",
    # repo = "https://github.zhaw.ch/hacs/AirfoilPolars",
    # sitename="AirfoilPolars.jl",
    modules=[AirfoilPolars],
    # format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Polar" => "polar.md",
        "Solve" => "solve.md",
        "IO" => "io.md",
    ]
)

deploydocs(;
    repo="https://github.zhaw.ch/hacs/AirfoilPolars.jl",
)