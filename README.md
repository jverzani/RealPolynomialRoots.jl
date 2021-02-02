# RealPolynomialRoots

[![Coverage](https://codecov.io/gh/jverzani/RealPolynomialRoots.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jverzani/RealPolynomialRoots.jl)


A package to find isolating intervals for the real roots of a square free polynomial.


Example:

```
julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1] # mignotte polynomial with two nearby roots

julia> ANewDsc(ps)
There were 3 isolating intervals found:
[3.0…, 9.5…]₅₃
[0.00787401589014…, 0.00787401713433…]₅₃
[0.00787401479283…, 0.00787401589014…]₅₃


julia> ps =[ # from https://discourse.julialang.org/t/root-isolation-of-real-rooted-integer-polynomials/51421
                      942438915208811912419937422298363203125
                   164182217245953398816894035758761902846875
                  4584900574568933770264468813466870772155175
                 48995332393110515074735075708247882042540865
                266674183150777010544241114017741621823207005
                852443280934837985352088128423887894438557515
               1738546146302892245736990844587487000484756535
               2381558158813900978436742173305983349418813145
               2262003889258248241081177038309445610985409335
               1516025051068561122302699855213604175083575145
                720810987764796866354279087485114863858085005
                241341213302325116160821160849326697595681275
                 55563538585119205063994483187179167406616375
                  8363912237256094118085070946083688310200625
                   740493466743082745510080711751444519503125
                    29215606371473169285018060091249259296875];

julia> ANewDsc(ps)
There were 15 isolating intervals found:
[-0.048706…, 0.0…]₁₇₃
[-0.10388…, -0.048706…]₁₇₃
[-0.1953…, -0.104…]₁₇₃
[-0.3691…, -0.1953…]₁₇₃
[-0.6953…, -0.3672…]₁₇₃
[-1.0…, -0.6953…]₁₇₃
[-1.344…, -1.0…]₁₇₃
[-1.688…, -1.344…]₁₇₃
[-2.062…, -1.688…]₁₇₃
[-2.469…, -2.062…]₁₇₃
[-2.844…, -2.469…]₁₇₃
[-3.25…, -2.844…]₁₇₃
[-3.594…, -3.25…]₁₇₃
[-3.797…, -3.594…]₁₇₃
[-4.0…, -3.797…]₁₇₃

```

The algorithm used is based on 

Computing Real Roots of Real Polynomials ... and now For Real!
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
arXiv:1605.00410; DOI:	10.1145/2930889.2930937


Performance could be **significantly** improved (polynomials of degree 10,000 or more are tractable with the algorithm, but this implementations only scales up to degree 100 or so). The `Hecke.jl` package has a much faster alternative.

