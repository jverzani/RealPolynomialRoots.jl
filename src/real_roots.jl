"""
    real_roots(p::Polynomials.StandardBasisPolynomial; square_free=true, [m], [M])

Find isolating intervals for the real roots of the univariate polynomial `p`.

* `p`: a polynomial in the standard basis
* `square_free`: The algorithm requires a square-free polynomial, which is assumed by default. If this argument is `false`, then a numerical GCD approach is used to identify the square free part (`p/gcd(p, p')`).
* `m`: an optional lower bound for the smallest possible real root
* `M`: an optional upper bound for the largest possible real root


Returns a "`State`" object which holds the isolating intervals in the component `Isol` and any unresolved intervals in the component `Unresolved`. Iterating over a `State` object, iterates over the identified isolating intervals.

Mostly follows the *basic* algorithm implemented in [`ANewDsc`](@ref) due to Alexander Kobel, Fabrice Rouillier, Michael Sagraloff [DOI](10.1145/2930889.2930937).


Example:

```jldoctest
julia> p = fromroots(Polynomial, [1,2,3]);

julia> real_roots(p)
There were 3 isolating intervals found:
[2.12…, 4.0…]₅₉
[1.19…, 2.12…]₅₉
[-0.0…, 1.19…]₅₉

julia> p = fromroots(Polynomial, [1,1,1,2,2,3]); # not square free

julia> real_roots(p, square_free=false)
There were 3 isolating intervals found:
[2.12…, 4.0…]₅₇
[1.19…, 2.12…]₅₇
[-0.0…, 1.19…]₅₇

julia> p = Polynomial([-1, 254, -16129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
Polynomial(-1 + 254*x - 16129*x^2 + x^15)

julia> real_roots(p)
There were 3 isolating intervals found:
[0.781…, 3.12…]₅₃
[0.00787401574803149570638…, 0.00787401574803149913348…]₂₁₂
[0.0078740157480314918219…, 0.0078740157480314957064…]₂₁₂

julia> filter(isreal, roots(p)) ## misses two with imaginary part ≈ 1e-10
1-element Array{Complex{Float64},1}:
 2.1057742291764914 + 0.0im
```

!!! note
    This implementation is 𝑶(n²); much slower than the `Hecke._roots` function provided through `arblib` in `Hecke.jl`, which itself says is not competitive with specialized algorithms, such as provided in the RS library of the paper authors.

!!! note
    The roots can be found from the intervals. See [`ANewDsc`](@ref) for an illustration
    using the `Roots` package to do so.
"""
function real_roots(p::Polynomials.StandardBasisPolynomial; square_free=true, kwargs...)
    if square_free
        return ANewDsc(coeffs(p); kwargs...)
    else
        u,v,w,_,_ = ngcd(p, derivative(p))
        return ANewDsc(coeffs(v); kwargs...)
    end
end

