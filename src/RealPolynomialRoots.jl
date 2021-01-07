module RealPolynomialRoots

using LinearAlgebra
using Random

import Roots
import Roots: find_zero

import ArbNumerics
import ArbNumerics: ArbFloat

include("realroots.jl")


# clean up
export Interval, State, anewdsc, nextstep!, upperbound, lowerbound, admissiblepoint
export descartesbound, arbL
export linearstep, poly_shift, zerotest, onetest, boundarytest, newtontest
   


end
