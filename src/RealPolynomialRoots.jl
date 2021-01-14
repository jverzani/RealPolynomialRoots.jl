module RealPolynomialRoots

using LinearAlgebra
using Random

import Roots
import Roots: find_zero

include("realroots.jl")

export ANewDsc, real_roots_sqfree

end
