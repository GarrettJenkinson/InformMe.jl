using Documenter, InformMe

makedocs(;
    modules=[InformMe],
    format=:html,
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
    latest = "master",
    target="build",
    deps = nothing,
    make = nothing,
    repo = "github.com/GarrettJenkinson/InformMe.jl",
    julia = "1.3",
    osname = "linux"
)
