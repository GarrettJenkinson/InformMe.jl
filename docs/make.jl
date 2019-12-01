using Documenter, InformMe

makedocs(;
    modules=[InformMe],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/GarrettJenkinson/InformMe.jl/blob/{commit}{path}#L{line}",
    sitename="InformMe.jl",
    authors="GarrettJenkinson",
    assets=String[],
)

deploydocs(;
    repo="github.com/GarrettJenkinson/InformMe.jl",
)
