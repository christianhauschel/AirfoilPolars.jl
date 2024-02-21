using Revise, Documenter, AirfoilPolars

makedocs(
    sitename="AirfoilPolars.jl",
    repo = "https://github.zhaw.ch/hacs/AirfoilPolars",
    # sitename="AirfoilPolars.jl",
    # modules=[AirfoilPolars],
    # format = Documenter.HTML(),
    pages = [
        # "polar.md",
        "Intro" => "index.md",
        "Polar" => "polar.md",
        "Solver" => "solve.md",
        # "Subsection" => [
        #     ...
        # ]
    ]
    
)