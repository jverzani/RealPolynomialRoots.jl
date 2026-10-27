# RealPolynomialRoots

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jverzani.github.io/RealPolynomialRoots.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/RealPolynomialRoots.jl/dev)
[![Build Status](https://github.com/jverzani/RealPolynomialRoots.jl/workflows/CI/badge.svg)](https://github.com/jverzani/RealPolynomialRoots.jl/actions)
[![codecov](https://codecov.io/gh/jverzani/RealPolynomialRoots.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jverzani/RealPolynomialRoots.jl)

A package to find isolating intervals for the real roots of a square free polynomial.

## Example

Polynomials are specified by their coefficients, ordered from the
constant term up. `ANewDsc` returns the isolating intervals:

```julia
julia> using RealPolynomialRoots

julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1]; # mignotte polynomial with two nearby roots

julia> st = ANewDsc(ps)
There were 3 isolating intervals found:
[3.5…, 7.75…]₅₃
[0.00787401546449…, 0.00787401667503…]₆₉
[0.00787401431353…, 0.00787401546449…]₆₉
```

The `refine_roots` function refines the intervals and returns approximate roots:

```julia
julia> refine_roots(st)
3-element Vector{BigFloat}:
 6.939437409621392124436713492447610272200680501712185816507667632045076114762801
 0.007874016089132754403608727898779727134193464194254785228308443233139694153163513
 0.007874015406930341157555003028161633376551552518768059431667490175426147456329789
```

## Algorithm

The algorithm used is based on:

*Computing Real Roots of Real Polynomials ... and now For Real!*
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
[arXiv](1605.00410) [DOI](https://doi.org/10.1145/2930889.2930937)

More detail on the algorithm is found in:

Computing real roots of real polynomials
Michael Sagraloff, Kurt Mehlhorn
[DOI](https://doi.org/10.1016/j.jsc.2015.03.004)
