module RealPolynomialRoots

using LinearAlgebra

include("bigfloat-utils.jl")
include("interval-arithmetic.jl")
include("descartes-bound.jl")
include("anewdsc.jl")
include("refine-roots.jl")

export ANewDsc, refine_roots

end
