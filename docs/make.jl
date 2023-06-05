using HdivBiotElasticityPaper
using Documenter

DocMeta.setdocmeta!(HdivBiotElasticityPaper, :DocTestSetup, :(using HdivBiotElasticityPaper); recursive=true)

makedocs(;
    modules=[HdivBiotElasticityPaper],
    authors="Santiago Badia <santiago.badia@monash.edu>, Martin Hornkjol <marhorn@math.uio.no>, Arbaz Khan <arbaz@ma.iitr.ac.in>, Kent-Andre Mardal <kent-and@simula.no>, Alberto F. Martin <alberto.f.martin@anu.edu.au>, Ricardo Ruiz-Baier <ricardo.ruizbaier@monash.edu>",
    repo="https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/blob/{commit}{path}#{line}",
    sitename="HdivBiotElasticityPaper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://amartinhuertas.github.io/HdivBiotElasticityPaper.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amartinhuertas/HdivBiotElasticityPaper.jl",
    devbranch="main",
)
