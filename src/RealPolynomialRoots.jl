module RealPolynomialRoots

using LinearAlgebra
using Random
using Polynomials
using Roots

using MutableArithmetics
const MA = MutableArithmetics

include("anewdsc.jl")
include("real_roots.jl")

export ANewDsc, real_roots

end
