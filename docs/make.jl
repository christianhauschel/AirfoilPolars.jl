using Documenter, AirfoilPolars

# push!(LOAD_PATH,"../src/")

makedocs(
    sitename="AirfoilPolars.jl",
    repo = "https://github.com/christianhauschel/AirfoilPolars.jl.git",
    # remotes=nothing,
    # modules=[AirfoilPolars],
    # format = Documenter.HTML(prettyurls=false),
    pages = [
        "Home" => "index.md",
        "Polar" => "polar.md",
        "Solve" => "solve.md",
        "IO" => "io.md",
    ]
)

deploydocs(;
    repo="github.com/christianhauschel/AirfoilPolars.jl.git",
)