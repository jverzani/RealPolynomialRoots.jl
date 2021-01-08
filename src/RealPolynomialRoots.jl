module RealPolynomialRoots

using LinearAlgebra
using Random

import Roots
import Roots: find_zero

import ArbNumerics
import ArbNumerics: ArbFloat

include("arb-utils.jl")
include("realroots.jl")

export anewdsc, real_roots_sqfree
# clean up
export Interval, State, anewdsc, nextstep!, _upperbound, _lowerbound, admissiblepoint
export descartesbound, arbL
export linearstep, poly_shift, zerotest, onetest, boundarytest, newtontest
   


end
