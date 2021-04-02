## Find isolating intervals for real roots of a polynomial
## --------------------------------------------------

## This implementation is too slow!
## While polys of degree 1000 or more are tractable, this is
## too slow to be usable for polynomials beyond degree 100 or so.
## The underlying implmentation here (not the main algorithm) does
## not scale as well as it should


const DEF_PRECISION = 53 


## Find precision of a big float
precision(x::Float64) = DEF_PRECISION
precision(x::BigFloat) = Base.MPFR.precision(x)

function poly_deriv(p)
    ntuple(i -> i*p[i+1], length(p)-1)
end

function poly_deriv!(p::Vector{T}) where {T}
    for i = 1:length(p)-1
        p[i] = i * p[i+1]
    end
    p[end] = zero(T)
    p
end


## The Taylor shift here is the most expensive operationm 𝑶(n²)
##
## https://arxiv.org/pdf/1605.00410.pdf has a better strategy
## of partial Taylor shifts, using just the nearby roots
##
##
## See Fast Approximate Polynomial Multipoint Evaluation and Applications (https://arxiv.org/pdf/1304.8069.pdf) for one possible speed up
## see  The Fundamental Theorem of Algebra in Terms of Computational Complexity
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.123.3313&rep=rep1&type=pdf
## for another
##
## compute [p(a), p'(a), 1/2p''(a). 1/3! p'''(a)...]
## as q(x) = p(x+a) = p(a) + p'(a)x + 1/2! p''(a)x^2 + 1/3! p'''(a) x^3 ....
## this uses O(n^2) Horner scheme
" `p(x + λ)`: Translate polynomial left by λ "
function poly_translate!(p::Vector{T}, λ::T=one(T)) where {T}
    n = length(p)
    dps = zeros(T, n)
    for i in n:-1:1
        for j in n:-1:2
            @inbounds dps[j] = muladd(dps[j], λ, dps[j-1])
        end
        @inbounds dps[1] = muladd(dps[1], λ, p[i])
    end

    copy!(p, dps)

    nothing

end

"Finds  `x^n p(1/x)` which is a reversal of coefficients "
function poly_reverse!(p)
    reverse!(p)
end

"scale x axis by λ"
function poly_scale!(p, λ)
    a = λ
    for i in 2:length(p)
        @inbounds p[i] *= a
        a *= λ
    end
end

"q = p(-x)"
function poly_flip!(p)
    for i in 2:2:length(p)
        @inbounds p[i] = -p[i]
    end
end

# shift polynomial so (a,b) -> (0, oo)
poly_shift(p, a, b) = poly_shift!(copy(p), a, b)
function poly_shift!(p::Vector{T}, a, b) where {T}
    poly_translate!(p, a)
    poly_scale!(p, b-a)
    poly_reverse!(p)
    poly_translate!(p, one(T))
    p
end

function descartesbound(p, a, b, L =  maximum(precision, (a,b)))
    T = BigFloat
    setprecision(L) do
        q = poly_shift(collect(T, p), T(a), T(b))
        descartescount(q, 1/T(2)^L)
    end
end
descartesbound(p) = descartescount(p, 0.0)

# Descartes bound on real roots based on the sign changes 
function descartescount(q, tol)
    cnt = -1
    s = zero(eltype(q))
    for qᵢ in q
        abs(qᵢ) < tol && return -1
        sᵢ = sign(qᵢ)
        if sᵢ != s
            s = sᵢ
            cnt += 1
        end
    end
    cnt
end

# m[ϵ] evenly spaced points in B(x,ϵ)
function multipoint(m, ϵ, n′)
    (m + i * ϵ/n′ for i in -n′:n′)
end

# This uses 𝒐(n^2) algorithm. This can be reduced
# The Alexander Kobel, Fabrice Rouillier, Michael Sagraloff paper suggests a randomization
# Schonhage has a method to compute in 𝒐(n ln(n)) time
function admissiblepoint(p, m, ϵ, n)

    T = BigFloat
    L = precision(m)

    n′ = sum(divrem(n,2))

    for i in 1:4
        tol = 1/big(2.0)^L
        out = setprecision(L) do
            tol = eps(BigFloat)
            mx, vx = zero(T), zero(T)
            for mᵢ in multipoint(T(m), T(ϵ), n)
                vᵢ = abs(evalpoly(mᵢ, p))
                if vᵢ > vx
                    vx, mx = vᵢ, mᵢ
                end
            end
            if vx > tol
                return T(mx)
            end
        end
        out != nothing && return out
        L = 2 * L
    end
    #@warn "no admissible point found"
    return m # pray?
end

## Hold an interval
struct Interval{T}
    a::T
    b::T
    N::Base.Ref{Int}
    Depth::Base.Ref{Int}
end

function Interval(a,b, n::Int=4, depth::Int=1)
    L = maximum(precision, (a,b))
    a′, b′ = setprecision(L) do
        big(a), big(b)
    end
    Interval(a′, b′, Ref(n), Ref(depth))
end

## An interval can be iterated over to return the end points
Base.length(I::Interval) = 2
function Base.iterate(I::Interval, state=nothing)
    if state==nothing
        return (I.a, 1)
    elseif state == 1
        return (I.b, 2)
    else
        return nothing
    end
end        

function Base.show(io::IO, I::Interval)
    a, b = I.a, I.b
    L = maximum(precision, (a,b))

    L′ = max(0, ceil(Int, -log2(b-a))) + 5
    ds = ceil(Int, -log10(b-a)) + 5
    
    sbs = ("₋","","","₀","₁","₂","₃","₄","₅","₆","₇","₈","₉")
    iob = IOBuffer()
    for i in string(L)
        print(iob, sbs[Int(i)-44])
    end
    ind = String(take!(iob))

    setprecision(L′) do
        print(io, "[$(round(a, digits=ds))…, $(round(b, digits=ds))…]$ind")
    end
end

function tₐ(pa)
    apa = abs(pa)
    L = precision(pa)
    apa > 1/big(2)^L && return ceil(Int, log2(apa))
    return 999999
end

τ(p) = max(1, log(norm(p, Inf)))
Mlog(x) = log(max(1, abs(x)))

# In https://people.mpi-inf.mpg.de/~msagralo/RealRootComputation.pdf 0-test
# the task of computing absolute L-bit approximations to p is not implemented
# as suggested, rather, we just evalute `descartesbound` using L bits of precision
# Similarly in onetest
function zerotest(p, a, b)

    L′ = maximum(precision, (a,b))
    pa, pb = setprecision(L′) do
        evalpoly(a, p), evalpoly(b, p)
    end
    ta, tb = tₐ(a), tₐ(b)

    n = length(p) + 1
    n′ = sum(divrem(n,2))

    L′′ = max(24, L′, max(1, -min(ta-1, tb-1) + 2(n+1) + 1)) 
    L = ceil(Int, n + τ(p) + n * log(max(1, abs(a))) + n * log(max(1, abs(b))) + L′′)

    mI = a + (b-a)/2
    bnda, bndb = descartesbound(p, a, mI, L), descartesbound(p, mI, b, L)
        
    iszero(bnda) && iszero(bndb) && return true
    return false
    
end


function onetest(p, a, b)

    L′ = maximum(precision, (a,b))
    ta, tb = tₐ(a), tₐ(b)

    n = length(p) - 1
    n′ = sum(divrem(n,2))

    ϵ = (b-a) / (4n)

    mstar = setprecision(L′) do
        admissiblepoint(p, a + (b-a)/2, ϵ, n′)
    end
    t = tₐ(mstar)
    
    L′′ = max(24, L′, max(1, -min(ta-1, tb-1, t-1) + 4(n+2))) 
    L = ceil(Int, n + τ(p) + n * log(max(1, abs(a))) + n * log(max(1, abs(b))) + L′′)

    bnda, bndb = descartesbound(p, a, mstar, L), descartesbound(p, mstar, b, L)
    
    iszero(bnda) && isone(bndb) && return true
    isone(bnda) && iszero(bndb) && return true

    return false

end

function newtontest(p, a, b, N)

    p, p′ = p, poly_deriv(p) 
    n = length(p) - 1
    n′ = sum(divrem(n, 2))
    ϵ = 1/(2^ceil(Int, 5 + log(n))) 

    L = maximum(precision, (a,b))

    for i in (1,2,3)
        for j in (i+1):3
            L′ = L
            for _ in 1:8
                val, J = setprecision(L′) do
                    m, w = a + (b-a)/2, b-a
                    ξᵢ, ξⱼ = admissiblepoint(p, a + i*w/4, ϵ * w, n′), admissiblepoint(p, a + j*w/4, ϵ * w, n′)
                    vᵢ, vⱼ = evalpoly(ξᵢ,p)/evalpoly(ξᵢ,p′), evalpoly(ξⱼ,p)/evalpoly(ξⱼ,p′)
                    
                    
                    if (abs(vᵢ) > w && abs(vⱼ) > w) ||
                        (abs(vᵢ - vⱼ) < w/(4n))
                        # (25) discard pair
                        return false, I
                    elseif abs(vᵢ) < 2*w && abs(vⱼ) < 2*w && abs(vᵢ-vⱼ) > w/(8n)
                        # (26)
                        k̄ = (ξⱼ - ξᵢ) / (vⱼ - vᵢ)
                        λᵢⱼ = ξᵢ - k̄ * vᵢ # k cluster is supposed
                        ## λ̃ is *supposed* to be a good approximation

                        a <= λᵢⱼ <= b  || return false, I
                        lᵢⱼ = floor(Int, Float64((λᵢⱼ - a) * (4 * N) / w))
                        
                        aᵢⱼ = a + max(0,  lᵢⱼ - 1) * w / (4N)
                        bᵢⱼ = a + min(4N, (lᵢⱼ + 2)) * w / (4N)

                        a⁺ᵢⱼ = (lᵢⱼ <= 1) ? a : admissiblepoint(p, aᵢⱼ, ϵ * w/N, n′)
                        (isnan(a⁺ᵢⱼ) || isinf(a⁺ᵢⱼ)) && return false, I
                        (!zerotest(p, a, a⁺ᵢⱼ)) && return false, I
                        b⁺ᵢⱼ = (bᵢⱼ == b) ? b : admissiblepoint(p, bᵢⱼ, ϵ * w/N, n′)
                        (isnan(b⁺ᵢⱼ) || isinf(b⁺ᵢⱼ)) && return false, I
                        (!zerotest(p, b⁺ᵢⱼ, b) ) && return false, I

                        return true, Interval(a⁺ᵢⱼ, b⁺ᵢⱼ)
                    end
                    return (nothing, I)
                end
                if val == nothing
                    L′ *= 2
                elseif val == false
                    # next i,j
                    continue
                else
                    return val, J
                end
            end
        end
    end
    return false, I
end


# retrun (nothing, I) -- no root in I
# return (true, J) can reduce interval
# return (false, I) can't trim, use linear
function boundarytest(p, a, b, N)


    n = length(p) - 1
    n′ = sum(divrem(n,2))
    ϵ = 1/2.0^(2 + ceil(Int, log(n)))

    L = maximum(precision, (a,b))

    setprecision(L) do

        w = b-a
        m = a + w/2
        
        mₗ, mᵣ = a + w/(2N), b - w/(2N)
        mₗ⁺, mᵣ⁺ = admissiblepoint(p, mₗ, ϵ*w/N, n′), admissiblepoint(p, mᵣ, ϵ*w/N, n′)
        a < mₗ⁺ <= mᵣ⁺ < b || return (false, I)

        zₗ, zᵣ = zerotest(p, mₗ⁺, b), zerotest(p, a, mᵣ⁺)
        
        zₗ && zᵣ && return (nothing, I) # no root
        zₗ && return (true, Interval(a, mₗ⁺))
        zᵣ && return (true, Interval(mᵣ⁺, b))
        (false, I)
    end

end

        
## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function upperbound(p)
    T = BigFloat
    #setprecision(DEF_PRECISION) do

    n = length(p) - 1
    L = 53 + n + ceil(Int, τ(p))
    setprecision(L) do
        p = p[findfirst(!iszero, p):end]
        descartesbound(p) == 0 && return zero(T)
        
        p = p[findfirst(!iszero, p):end]
        
        q, d = p/p[end], length(p)-1
        
        d == 0 && error("degree 0 is a constant")
        d == 1 && abs(q[1])
        
        
        a1 = abs(q[d])
        B = maximum([abs(q[i]) for i in 1:(d-1)])
        
        a,b,c = 1, -(1+a1), a1-B
        out = (-b + sqrt(b^2 - 4a*c))/2
        
        out
    end
end

function lowerbound(p)
    q = p[findfirst(!iszero, p):end]
    poly_flip!(q)

    n = length(p) - 1
    L = 53 + n + ceil(Int, τ(p))
    setprecision(L) do
        -upperbound(q)
    end
end

# hold the state
struct State{N,T}
    Isol::Vector{Interval{BigFloat}}                        # DesBound == 1
    Unresolved::Vector{Interval{BigFloat}}
    p::NTuple{N,T}
end

(st::State)(x) = evalpoly(x, st.p)

## iterate over Isol
Base.length(st::State) = length(st.Isol)
function Base.iterate(st::State, state=nothing)
    if state==nothing
        return iterate(st.Isol)
    else
        return iterate(st.Isol, state)
    end
end        

function Base.show(io::IO, st::State)
    if !isempty(st.Unresolved)
        println(io, "There are unresolved intervals:")
        for I in st.Unresolved
            println(io, I)
        end
        println(io, "")
    end

    Is = st.Isol
    n = length(Is)

    if n == 0
        println(io, "No isolating intervals found")
        return nothing
    elseif n == 1
        println(io, "There was 1 isolating interval found:")
    else
        println(io, "There were $(length(Is)) isolating intervals found:")
    end
    
    for I in Is
        println(io, I)
    end
    
end

"""
    ANewDsc(p; [m=lowerbound(p)], [M=upperbound(p)])

A method to find isolating intervals for the real roots of the polynomial specified by `p`.

* `p`: the polynomial coefficients, `[a₀, a₁, …, aₙ]`, of a **square-free** polynomial.
* `m`: a lower bound for the smallest possible real root
* `M`: an upper bound for the largest possible real root

Returns a `State` instance which has components:

* `Isol` holding the isolating intervals. Iteration over a `State` object will iterate over `Isol`.
* `Unresolved` holding any unresolved intervals. The show method alerts the presence of any such intervals.


Examples:

```jldoctest
julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1] # mignotte polynomial with two nearby roots

julia> st = ANewDsc(ps)
There were 3 isolating intervals found:
[3.0…, 9.5…]₅₃
[0.00787401589014…, 0.00787401713433…]₅₃
[0.00787401479283…, 0.00787401589014…]₅₃

julia> ps = [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]; # πᵢ₌₁¹⁰ (x-i)

julia> st = ANewDsc(ps)
There were 10 isolating intervals found:
[9.5…, 11.5…]₈₀
[8.75…, 9.5…]₈₀
[7.75…, 8.5…]₈₀
[6.75…, 7.75…]₈₀
[5.25…, 6.75…]₈₀
[4.5…, 5.25…]₈₀
[3.56…, 4.5…]₈₀
[2.62…, 3.5…]₈₀
[1.69…, 2.56…]₈₀
[-0.0…, 1.69…]₈₀

julia> ps =[ # from https://discourse.julialang.org/t/root-isolation-of-real-rooted-integer-polynomials/51421/1
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

julia> ANewDsc(ps, m=-4.0, M=0.0)
There were 15 isolating intervals found:
[-0.053772…, 0.0…]₅₃
[-0.11475…, -0.053772…]₅₃
[-0.2451…, -0.1147…]₅₃
[...]
```

To find the zeros, the `Roots` package may be used. For example:

```
julia> using Roots

julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1];

julia> st = ANewDsc(ps)
There were 3 isolating intervals found:
[3.0…, 9.5…]₅₃
[0.00787401589014…, 0.00787401713433…]₅₃
[0.00787401479283…, 0.00787401589014…]₅₃

julia> [find_zero(st, I, Roots.BisectionExact()) for I ∈ st]
3-element Array{BigFloat,1}:
 6.93943740962139212443671349244761027220068050171218581650766763204507611476287
 0.00787401608913275440360872789877972713419346419425478522830844323313969358762523
 0.007874015406930341157555003028161633376551552518768059431667490175426147375209448
```

The default bracketing method for `BigFloat` can have issues with some problems, so we use `BisectionExact` above. Alternatively, `Roots.bisection(st, I...)` may be used.

Comparing to some alternatives, we have, the functionality from
Hecke.jl (`Hecke._roots`) is much better. However, compared to others this
could be seen as useful:

```
julia> x = variable(Polynomial);

julia> p = -1 + 254*x - 16129*x^2 + x^15;

julia> @time real_roots(p) # ANewDsc(coeffs(p))
  0.091078 seconds (2.41 M allocations: 106.406 MiB, 27.02% gc time)
There were 3 isolating intervals found:
[0.781…, 3.12…]₅₃
[0.00787401574803149570638…, 0.00787401574803149913348…]₂₁₂
[0.0078740157480314918219…, 0.0078740157480314957064…]₂₁₂

julia> filter(isreal, roots(p)) # much faster, but misses two roots with imaginary part ~ 1e-10
1-element Array{Complex{Float64},1}:
 2.1057742291764914 + 0.0im

julia> filter(isreal, AMRVW.roots(Float64.(coeffs(p)))) # using AMRVW. Similarly misses two roots
1-element Array{Complex{Float64},1}:
 2.1057742291764407 + 0.0im

julia> filter(isreal, PolynomialRoots.roots(coeffs(p))) # using PolynomialRoots. Misses 2.105?
2-element Array{Complex{Float64},1}:
   0.0078740158234482 + 0.0im
 0.007874015672614794 + 0.0im

julia> IntervalRootFinding.roots(f, IntervalArithmetic.Interval(0.0, 5.0)) # using IntervalRootFinding, IntervalArithmetic
8-element Array{Root{IntervalArithmetic.Interval{Float64}},1}:
 Root([0.00787395, 0.00787397], :unknown)
 Root([0.00787303, 0.00787305], :unknown)
 Root([0.00787418, 0.00787422], :unknown)
 Root([0.00787403, 0.00787409], :unknown)
 Root([2.10577, 2.10578], :unique)
 Root([0.00787383, 0.0078739], :unknown)
 Root([0.00787396, 0.00787404], :unknown)
 Root([0.00787357, 0.00787361], :unknown)

julia> @vars x # using SymPy

julia> @time  rts = sympy.real_roots(p(x)); # correctly identifies 3. 
  0.003896 seconds (518 allocations: 13.359 KiB)

julia> sympy.N(rts[2]) # takes a long time! (162s)
0.00787401574803150
```

The algorithm used is a partial implementation of

Computing Real Roots of Real Polynomials ... and now For Real!
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
arXiv:1605.00410; DOI:	10.1145/2930889.2930937

The algorithm relies on Descartes' rule of signs, which gives a bound
on the number of positive real roots of a polynomial, `p(x)`. By
considering the polynomial `p((ax+b)/(x+1))` (a mapping of `[a,b]` to
`[0,∞)` a bound on the number of roots in `[a,b]` can be found. A
simple bisection algorithm for a square-free polynomial would be to
start with an interval large enough to contain all the real roots,
then subdivide at the midpoint throwing out sub intervals which are
found to have no root; saving intervals known to have 1 root, and
repeating the subdivision otherwise. Issues with this are the need for
many subdivisions when clusters of roots are present and the numeric
issues that arise in computing the mapping.

The work of Kobel, Rouillier, and Sagraloff improves this by
introducing a Newton test for rapidly decreasing the size of an
interval when possible, and the ability to use finite precision
arithmetic, instead of exact arithmetic or interval arithmetic, to
compute the Descartes' bound, in addition to other algorithmic
improvements (not implemented here).

!!! note
    This implementation is 𝑶(n²); much slower than the `Hecke._roots`
    function provided through `arblib` in `Hecke.jl`, which itself
    says is not competitive with more specialized algorithms, such as
    provided in the RS library of the paper authors.

!!! note
    A square free polynomial can be found through `p/gcd(p, p')`,
    though in practice this calculation is numerically unstable.


"""
function ANewDsc(p; m=lowerbound(p), M=upperbound(p), max_depth=96)

    DEBUG = false
    
    n = length(p) - 1
    n′ = sum(divrem(n,2))

    L = maximum(precision∘float, (m, M))
    m′, M′ = setprecision(L) do
        big(m), big(M)
    end

    
    I = Interval(m′, M′)
    I.N[], I.Depth[] = 4, 1
    Internal = [I]
    Isol = typeof(I)[]
    Unresolved = typeof(I)[]

    while !isempty(Internal)
        I = pop!(Internal)
        a,b = I
        
        # are we resolved sufficiently
        if zerotest(p, a, b)
            continue
        end

        if onetest(p, a, b)
            push!(Isol, I)
            continue
        end

        
        ## pump the brakes, if needed
        if I.Depth[] > max_depth
            DEBUG && @show :max_depth
            push!(Unresolved, I)
            continue
        end


        L = maximum(precision, (a, b))        
        
        pa, pb = setprecision(L) do
            evalpoly(a, p), evalpoly(b, p)
        end

        if iszero(pa) || iszero(pb)
            DEBUG && @show :early_zero
            push!(Unresolved, I)
            continue
        end
        
        # shrink or divide
        N, depth′ = I.N[], I.Depth[] + 1

        val, J = newtontest(p, a, b, N)
        if val
            λ = Float64((I.b-I.a)/(J.b - J.a))
            DEBUG && @show :newton, λ
            J.N[], J.Depth[] = N^2, depth′ 
            push!(Internal, J)
            continue
        end
        
        val, J = boundarytest(p, a, b, N)
        if val == nothing
            # no root
            continue
        elseif val
            DEBUG && @show :boundary
            J.N[], J.Depth[] = N, depth′             
            push!(Internal, J)
            continue
        end


        setprecision(L) do
            w = b-a
            m = a + (b-a)/2
            c = admissiblepoint(p, m, w/32, n′)
            N′ = max(4, ceil(Int, sqrt(N)))
            DEBUG && @show :linear
            push!(Internal, Interval(a,c, N′, depth′))
            push!(Internal, Interval(c,b, N′, depth′))
        end
    end
    
    return State(Isol, Unresolved, NTuple{n+1,eltype(p)}(p))
    

end
