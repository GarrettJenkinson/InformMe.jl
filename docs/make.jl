using Documenter, InformMe

makedocs(;
    modules=[InformMe],
    format=Documenter.HTML(),
    source="src",
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/GarrettJenkinson/InformMe.jl/blob/{commit}{path}#L{line}",
    sitename="InformMe.jl",
    authors="Garrett Jenkinson"
)

deploydocs(;
    branch = "gh-pages",
    target="build",
    repo = "github.com/GarrettJenkinson/InformMe.jl"
)
