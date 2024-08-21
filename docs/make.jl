using ThreeBodyDecays
using Documenter
using Literate

DocMeta.setdocmeta!(
    ThreeBodyDecays,
    :DocTestSetup,
    :(using ThreeBodyDecays);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end

gen_content_dir = joinpath(@__DIR__, "src")
name = "10-tutorial"
tutorial_src = joinpath(@__DIR__, "src", "$(name).jl")
Literate.markdown(tutorial_src, gen_content_dir; name, documenter=true, credit=true, postprocess=fix_literate_output)

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
