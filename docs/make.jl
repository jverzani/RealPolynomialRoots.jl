using Documenter
using RealPolynomialRoots


ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"


DocMeta.setdocmeta!(RealPolynomialRoots, :DocTestSetup, :(using RealPolynomialRoots); recursive = true)

makedocs(
    sitename = "RealPolynomialRoots",
    format = Documenter.HTML(ansicolor=true),
    modules = [RealPolynomialRoots]
)

deploydocs(
    repo = "github.com/jverzani/RealPolynomialRoots.jl"
)
