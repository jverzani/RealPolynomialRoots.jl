"""
    real_roots(p::Polynomials.StandardBasisPolynomial; square_free=true, [m], [M])

Find the real roots of the univariate polynomial `p`.

* `p`: a polynomial in the standard basis
* `square_free`: The algorithm requires a square-free polynomial, which is assumed by default. If this argument is `false`, then a numerical GCD approach is used to identify the square free part (`p/gcd(p, p')`).
* `m`: an optional lower bound for the smallest possible real root
* `M`: an optional upper bound for the largest possible real root

Returns a container of identified real roots.

Use the [`ANewDsc`](@ref) function, 
which mostly follows the basic implementation of an  algorithm due to Alexander Kobel, Fabrice Rouillier, Michael Sagraloff [DOI](10.1145/2930889.2930937). The implementation here is *much* slower than others.


Example:

```jldoctest
julia> p = fromroots(Polynomial, [1,2,3]);

julia> real_roots(p)
3-element Vector{Float64}:
 3.0
 2.0
 1.0

julia> p = fromroots(Polynomial, [1,1,1,2,2,3]); # not square free

julia> real_roots(p, square_free=false)
3-element Vector{Float64}:
 3.0000000000041767
 2.000000000000601
 0.9999999999999366

julia> p = Polynomial([-1, 254, -16129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
Polynomial(-1 + 254*x - 16129*x^2 + x^15)

julia> real_roots(p)
3-element Vector{AbstractFloat}:
 2.1057742291764834
 0.007874015748031497374190409031015351762713667398747530835438348975
 0.0078740157480314949573666522672607058737133311592752406327723222706

julia> filter(isreal, roots(p)) ## misses two with imaginary part ≈ 1e-10
1-element Array{Complex{Float64},1}:
 2.1057742291764914 + 0.0im
```

!!! note
    This implementation is 𝑶(n²); much slower than the `Hecke._roots` function provided through `arblib` in `Hecke.jl`, which itself says is not competitive with specialized algorithms, such as provided in the RS library of the paper authors.

"""
function real_roots(p::P; square_free=true, kwargs...) where {T, P <: Polynomials.StandardBasisPolynomial{T}}
    if square_free
        st = ANewDsc(p; kwargs...)
    else
        u,v,w,_,_ = Polynomials.ngcd(p, derivative(p))
        st = ANewDsc(v; kwargs...)
    end

    if !isempty(st.Unresolved)
        @warn "there are unresolved intervals"
    end

    map(st) do I
        prec = precision(I.a)
        if prec == 53
            find_zero(st, Float64.(I))
        else
            setprecision(prec) do
                find_zero(st, I, Roots.BisectionExact())
            end
        end
    end
    
end

