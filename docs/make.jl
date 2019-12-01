using Documenter, InformMe

makedocs(;
    modules=[InformMe],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/GarrettJenkinson/InformMe.jl/blob/{commit}{path}#L{line}",
    sitename="InformMe.jl",
    authors="GarrettJenkinson"
)

deploydocs(;
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/GarrettJenkinson/InformMe.jl",
    julia = "1.3",
    osname = "linux"
)
