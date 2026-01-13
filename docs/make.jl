using Jspec
using Documenter

DocMeta.setdocmeta!(Jspec, :DocTestSetup, :(using Jspec); recursive=true)

makedocs(;
         modules=[Jspec],
    authors="Stefano Covino <stefano.covino@inaf.it> and contributors",
    sitename="Jspec.jl",
    format=Documenter.HTML(;
                           canonical="https://stefanocovino.github.io/Jspec.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
           repo="github.com/stefanocovino/Jspec.jl",
    devbranch="main",
)
