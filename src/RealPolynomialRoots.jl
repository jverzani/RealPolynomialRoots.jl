module RealPolynomialRoots

using LinearAlgebra
using Random

import Roots
import Roots: find_zero

using ArbNumerics


include("arb-utils.jl")
include("realroots.jl")

export anewdsc, real_roots_sqfree

end
