using ThreeBodyDecays
using Documenter

DocMeta.setdocmeta!(
    ThreeBodyDecays,
    :DocTestSetup,
    :(using ThreeBodyDecays);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
    modules = [ThreeBodyDecays],
    authors = "Misha Mikhasenko <mikhail.mikhasenko@cern.ch> and contributors",
    repo = "https://github.com/mmikhasenko/ThreeBodyDecays.jl/blob/{commit}{path}#{line}",
    sitename = "ThreeBodyDecays.jl",
    format = Documenter.HTML(;
        canonical = "https://mmikhasenko.github.io/ThreeBodyDecays.jl",
        repolink = "https://github.com/mmikhasenko/ThreeBodyDecays.jl",
    ),
    pages = [
        "index.md"
        [
            file for file in readdir(joinpath(@__DIR__, "src")) if
            file != "index.md" && splitext(file)[2] == ".md"
        ]
    ],
)

deploydocs(; repo = "github.com/mmikhasenko/ThreeBodyDecays.jl")
