# RealPolynomialRoots

[(source) RealPolynomialRoots](https://github.com/jverzani/RealPolynomialRoots.jl/)


A package to find isolating intervals for the real roots of a square free polynomial.


Example:

```@repl RPR
using RealPolynomialRoots
ps = [-1, 254, -16129, 0, 0, 0, 0, 1] # mignotte polynomial with two nearby roots
ANewDsc(ps)
```


```@repl RPR
ps =[ # from https://discourse.julialang.org/t/root-isolation-of-real-rooted-integer-polynomials/51421
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
ANewDsc(ps)
```

The `refine_roots` method refines the isolating intervals down to
diameter `1/2^L` and then takes the midpoint:

```@repl RPR
ps = [-1, 254, -16129, 0, 0, 0, 0, 1];
(refine_roots ∘ ANewDsc)(ps)
```

The algorithm used is based on

*Computing Real Roots of Real Polynomials ... and now For Real!*
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
([arXiv](https://arXiv.org/1605.00410), [DOI](https://doi.org/10.1145/2930889.2930937)).

More detail on the algorithm is found in

*Computing real roots of real polynomials*
Michael Sagraloff, Kurt Mehlhorn
[DOI](https://doi.org/10.1016/j.jsc.2015.03.004)



The performance of this implementation could be **significantly**
improved upon: polynomials of degree 10,000 or more are tractable with
the algorithm, but this implementation gets pretty sluggish on random
polynomials of degree 250.


## Reference

```@autodocs
Modules = [RealPolynomialRoots]
```
