using Documenter
using MaterialDocs
using RealPolynomialRoots


ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"


DocMeta.setdocmeta!(RealPolynomialRoots, :DocTestSetup, :(using RealPolynomialRoots); recursive = true)

makedocs(
    sitename = "RealPolynomialRoots",
    format = Material3(ansicolor=true, dark_mode=:toggle),
    modules = [RealPolynomialRoots]
)

deploydocs(
    repo = "github.com/jverzani/RealPolynomialRoots.jl"
)
