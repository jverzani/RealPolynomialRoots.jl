## ----
Container{T} = Union{AbstractVector{T}, NTuple{N,T}} where {N}

## ---- Interval

## Hold an interval
struct 𝑰
    a::BigFloat
    b::BigFloat
    N::Base.Ref{Int}
    Depth::Base.Ref{Int}
end

function 𝑰(a,b, n::Int=4, depth::Int=1)
    L = maximum(precision, (a,b))
    a′, b′ = setprecision(L) do
        big(a), big(b)
    end
    𝑰(a′, b′, Ref(n), Ref(depth))
end

## An interval can be iterated over to return the end points
Base.length(I::𝑰) = 2
function Base.iterate(I::𝑰, state=nothing)
    isnothing(state) && return (I.a, 1)
    isone(state) && return (I.b, 2)
    return nothing
end        

function Base.show(io::IO, I::𝑰)
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

## ---- bigfloat speedups
## These use MutableArithmetics.jl
using MutableArithmetics
const MA = MutableArithmetics

# a,b -> b,a using temporary c
function swap!(c, a,b)
    MA.zero!(c)
    MA.add!(c, b)
    MA.zero!(b)
    MA.add!(b, a)
    MA.zero!(a)
    MA.add!(a, c)
end

# a <- a*b+c
function muladd!(a::T,b::T,c::T) where {T <: BigFloat}
    MA.mul!(a,b)
    MA.add!(a,c)
end

# allocates less
function _evalpoly(x, ps::Container{T}) where {T}
    N = length(ps)
    ex = one(T)
    for i in N-1:-1:1
        muladd!(ex, x, ps[i])
    end
    ex
end

## --- basic transformations

# p -> p(-x)
function Base.reverse!(p::NTuple{N,T}) where {N,T}
    n = N ÷ 2
    c = one(T)
    for i ∈ 1:n
        j = N + 1 - i
        pᵢ, pⱼ = p[i], p[j]
        swap!(c, pᵢ,pⱼ)
    end
    p
end

# p -> p(λx)
function scale!(p::Container{T}, λ::T) where {T}
    a = one(T)
    n = length(p)
    for i ∈ 2:n
        pᵢ = p[i]
        MA.mul!(a, λ)
        MA.mul!(pᵢ, a)
    end
    return nothing
end

# p -> p(x + λ)  This is 𝑶(n²) Horner scheme and could be faster
function taylor_shift!(ps::Container{T}, λ=one(T)) where {T}
    N = length(ps)
    dps = deepcopy(ps)
    MA.zero!.(ps)
    @inbounds for i in N:-1:1
        for j in N:-1:2
            muladd!(ps[j], λ, ps[j-1])
        end
        muladd!(ps[1], λ, dps[i])
    end
    return nothing
end

# p -> p((ax + b)/(x+b))
function mobius_transform!(p::Container{T}, a, b) where {T}
    taylor_shift!(p, a)
    scale!(p, (b-a))
    reverse!(p)
    taylor_shift!(p, one(T))
    nothing
end      


## -----

const DEF_PRECISION = 53 

## Find precision of a big float
precision(x::Float64) = DEF_PRECISION
precision(x::BigFloat) = Base.MPFR.precision(x)

## -----

function descartesbound(p, a, b, L =  maximum(precision, (a,b)))
    T = BigFloat
    setprecision(L) do
        q = T.(deepcopy.(p))
        mobius_transform!(q, T(a), T(b))        
        u = descartescount(q, 1/T(2)^L)
        u
    end
end
descartesbound(p::Container{T}) where {T}  = descartescount(p, zero(T))

# Descartes bound on real roots based on the sign changes 
function descartescount(q::Container{T}, tol::T) where {T}
    n = length(q)
    q₀ = q[1]

    tol_ = -tol
    flag = q₀ < tol ? true : false
    flag && q₀ > tol_ && return -1

    cnt = 0

    for i ∈ 2:n
        qᵢ = q[i]
        flag′ = qᵢ < tol ? true : false
        flag′ && qᵢ > tol_ && return -1
        if flag′ != flag
            cnt += 1
            flag = flag′
        end
    end
    cnt
end

## ----

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
                vᵢ = abs(_evalpoly(mᵢ, p))
                if vᵢ > vx
                    vx, mx = vᵢ, mᵢ
                end
            end
            if vx > tol
                return T(mx)
            end
            return nothing
        end
        out != nothing && return out
        L = 2 * L
    end
    #@warn "no admissible point found"
    return T(m) # pray?
end

## ---- Tests (onetest, zerotest, newtontest, boundarytest)

function tₐ(pa)
    apa = abs(pa)
    L = precision(pa)
    apa > 1/big(2)^L && return ceil(Int, log2(apa))
    return 999999
end

τ(p::Container{T}) where {T} = max(one(T), log(norm(p, Inf)))
Mlog(x) = log(max(1, abs(x)))

# In https://people.mpi-inf.mpg.de/~msagralo/RealRootComputation.pdf 0-test
# the task of computing absolute L-bit approximations to p is not implemented
# as suggested, rather, we just evalute `descartesbound` using L bits of precision
# Similarly in onetest
function zerotest(p, a, b)
    a >= b && return true
    L′ = maximum(precision, (a,b))
    pa, pb = setprecision(L′) do
        _evalpoly(a, p), _evalpoly(b, p)
        #p(a), p(b)
    end
    ta, tb = tₐ(a), tₐ(b)

    n = length(p) - 1
    n′ = sum(divrem(n,2))

    L′′ = max(24, L′, max(1, -min(ta-1, tb-1) + 2(n+1) + 1)) 
    L = ceil(Int, n + τ(p) + n * log(max(1, abs(a))) + n * log(max(1, abs(b))) + L′′)

    mI = a + (b-a)/2
    bnda, bndb = descartesbound(p, a, mI, L), descartesbound(p, mI, b, L)

    iszero(bnda) && iszero(bndb) && return true
    return false
    
end

function onetest(p, a, b)
    T = BigFloat
    a >= b && return false
    L′ = maximum(precision, (a,b))
    ta, tb = tₐ(a), tₐ(b)

    n = length(p) - 1
    n′ = sum(divrem(n,2))

    ϵ = (b-a) / (4n)

    mstar = setprecision(L′) do
        admissiblepoint(p, a + (b-a)/2, ϵ, n′)
    end
    t::T = tₐ(mstar)
    
    L′′ = max(24, L′, max(1, -min(ta-1, tb-1, t-1) + 4(n+2))) 
    L = ceil(Int, n + τ(p) + n * log(max(1, abs(a))) + n * log(max(1, abs(b))) + L′′)

    bnda, bndb = descartesbound(p, a, mstar, L), descartesbound(p, mstar, b, L)
    
    iszero(bnda) && isone(bndb) && return true
    isone(bnda) && iszero(bndb) && return true

    return false

end

function newtontest(p, p′, a, b, N)
    n = length(p) - 1
    n′ = sum(divrem(n, 2))
    ϵ = 1/(2^ceil(Int, 5 + log(n))) 
    I = 𝑰(a,b,N)
    L = maximum(precision, (a,b))

    for i in (1,2,3)
        for j in (i+1):3
            L′ = L
            for _ in 1:8
                val, J = setprecision(L′) do
                    m, w = a + (b-a)/2, b-a
                    ξᵢ, ξⱼ = admissiblepoint(p, a + i*w/4, ϵ * w, n′), admissiblepoint(p, a + j*w/4, ϵ * w, n′)
                    vᵢ,vⱼ = _evalpoly(ξᵢ,p)/_evalpoly(ξᵢ,p′), _evalpoly(ξⱼ,p)/_evalpoly(ξⱼ,p′)
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
                        !zerotest(p, a, a⁺ᵢⱼ) && return false, I
                        b⁺ᵢⱼ = (bᵢⱼ == b) ? b : admissiblepoint(p, bᵢⱼ, ϵ * w/N, n′)
                        (isnan(b⁺ᵢⱼ) || isinf(b⁺ᵢⱼ)) && return false, I
                        !zerotest(p, b⁺ᵢⱼ, b) && return false, I
                        return true, 𝑰(a⁺ᵢⱼ, b⁺ᵢⱼ)
                    end
                    #return (nothing, I)
                    return (false, I)                    
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


# retrun (0, I) -- no root in I
# return (1, J) can reduce interval
# return (-1, I) can't trim, use linear
function boundarytest(p, a, b, N)


    n = length(p) - 1
    n′ = sum(divrem(n,2))
    ϵ = 1/2.0^(2 + ceil(Int, log(n)))
    I = 𝑰(a,b,N)
    L = maximum(precision, (a,b))

    setprecision(L) do

        w = b-a
        m = a + w/2
        
        mₗ, mᵣ = a + w/(2N), b - w/(2N)
        mₗ⁺, mᵣ⁺ = admissiblepoint(p, mₗ, ϵ*w/N, n′), admissiblepoint(p, mᵣ, ϵ*w/N, n′)
        a < mₗ⁺ <= mᵣ⁺ < b || return (-1, I)
        
        zₗ, zᵣ = zerotest(p, mₗ⁺, b), zerotest(p, a, mᵣ⁺)
        
        zₗ && zᵣ && return (0, I) # no root        
        zₗ && return (1, 𝑰(a, mₗ⁺))
        zᵣ && return (1, 𝑰(mᵣ⁺, b))
        (-1, I)
    end

end

## ----

## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function upperbound(p)
    T = BigFloat
    n = length(p) - 1
    
    L = 53 + n + ceil(Int, τ(p))
    setprecision(L) do

        descartesbound(p) == 0 && return zero(T)
        q, d = p./p[end], n
        
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
    q = deepcopy(p)
    scale!(q, -one(eltype(q)))
    d = length(q) - 1
    d <= 0 && error("constant")

    L = 53 + d + ceil(Int, τ(p))
    setprecision(L) do
        -upperbound(q)
    end
end

## ----

# hold the state
struct State{N,T}
    Isol::Vector{𝑰}                        # DesBound == 1
    Unresolved::Vector{𝑰}
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

## ----

"""
    ANewDsc(p; [m=lowerbound(p)], [M=upperbound(p)])

A method to find isolating intervals for the real roots of a square-free polynomial specified by `p`.

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
[4.25…, 17.0…]₅₃
[0.00787401572688…, 0.00787401786147…]₅₃
[0.00787401345315…, 0.00787401572688…]₅₃

julia> ps = [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]; # πᵢ₌₁¹⁰ (x-i)

julia> st = ANewDsc(ps)
There were 10 isolating intervals found:
[9.75…, 10.2…]₅₃
[9.0…, 9.75…]₅₃
[7.75…, 9.0…]₅₃
[6.0…, 7.75…]₅₃
[5.38…, 6.12…]₅₃
[4.5…, 5.38…]₅₃
[3.12…, 4.5…]₅₃
[2.5…, 3.19…]₅₃
[1.81…, 2.5…]₅₃
[-3.38…, 1.81…]₅₃

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

julia> ANewDsc(ps)
There were 15 isolating intervals found:
[-0.053406…, 0.0…]₁₇₃
[-0.11389…, -0.053406…]₁₇₃
[-0.2383…, -0.1138…]₁₇₃
[-0.4492…, -0.2383…]₁₇₃
[-0.6602…, -0.4492…]₁₇₃
[...]
```

Comparing to some alternatives, we have that the functionality from
Hecke.jl (`Hecke._roots`) is **much** better. 

However, compared to other alternatives this could be seen as useful:

```
julia> x = variable(Polynomial);

julia> p = -1 + 254*x - 16129*x^2 + x^15;

julia> @time real_roots(p)
  0.114760 seconds (664.03 k allocations: 35.719 MiB, 14.09% gc time)
3-element Vector{AbstractFloat}:
 2.1057742291764834
 0.007874015748031497374190409031015351762713667398747530835438348975
 0.0078740157480314949573666522672607058737133311592752406327723222706

julia> filter(isreal, roots(p)) # much faster, but misses two roots with imaginary part ~ 1e-10
1-element Array{Complex{Float64},1}:
 2.1057742291764914 + 0.0im

julia> filter(isreal, AMRVW.roots(Float64.(coeffs(p)))) # using AMRVW. Similarly misses two roots
1-element Array{Complex{Float64},1}:
 2.1057742291764407 + 0.0im

julia> filter(isreal, AMRVW.roots(BigFloat.(coeffs(p)))) # this works here
3-element Vector{Complex{BigFloat}}:
 0.007874015748031494751793842937491860399146218747427882112208862187156603046408902 + 0.0im
 0.007874015748031497374190409031015351762713667398750832139185145345899098824178322 + 0.0im
    2.105774229176482954331883383107195997983004314462374928263620342390986189650864 + 0.0im

julia> filter(isreal, PolynomialRoots.roots(coeffs(p))) # using PolynomialRoots. Misses 2.105?
2-element Array{Complex{Float64},1}:
   0.0078740158234482 + 0.0im
 0.007874015672614794 + 0.0im

julia> IntervalRootFinding.roots(x->p(x), IntervalArithmetic.Interval(0.0, 5.0)) # using IntervalRootFinding, IntervalArithmetic
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
    A square free polynomial can be found through `p/gcd(p, p')`,
    though in practice this calculation is numerically unstable.

!!! note This implementation is **much** slower than the
    `Hecke._roots` function provided through `arblib` in `Hecke.jl`,
    which itself says is not competitive with more specialized
    algorithms, such as provided in the RS library of the paper
    authors. There are several reasons: The `mobius_transform!`
    function is 𝑶(n²), and could be 𝑶(n⋅log(n)) with more effort; the
    polynomial evaluation in `admissiblepoint` could, similarly, be
    made more efficient; despite using `MutableArithmetics.jl` to
    reduce allocations with the `BigFloat` type, there are still *far*
    too many allocations; the significant engineering speedups
    suggested by Kobel, Rouillier, and Sagraloff are not implemented;
    etc. 

    This implementation also is not as careful with floating point
    considerations as needed and detailed in the paper. The broken
    test with separated roots ~ 1e-100 illustrates this.



"""
function ANewDsc(q; m=lowerbound(q), M=upperbound(q), max_depth=96)
    
    p = BigFloat.(q)
    n = length(p) - 1
    n′ = sum(divrem(n,2))
    p′ = [i*p[i+1] for i ∈ 1:n]

    L = maximum(precision∘float, (m, M))
    m′, M′ = setprecision(L) do
        big(m), big(M)
    end

    
    I = 𝑰(m′, M′)
    I.N[], I.Depth[] = 4, 1
    Internal = [I]
    Isol = typeof(I)[]
    Unresolved = typeof(I)[]
    ctr = 0
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
            push!(Unresolved, I)
            continue
        end


        L = maximum(precision, (a, b))        
        
        pa, pb = setprecision(L) do
            _evalpoly(a,p), _evalpoly(b,p)
        end

        if iszero(pa) || iszero(pb)
            push!(Unresolved, I)
            continue
        end
        
        # shrink or divide
        N, depth′ = I.N[], I.Depth[] + 1
        ctr += 1

        val, J = newtontest(p, p′, a, b, N)
        if val
            λ = Float64((I.b-I.a)/(J.b - J.a))
            J.N[], J.Depth[] = N^2, depth′ 
            push!(Internal, J)
            continue
        end
        
        val, J = boundarytest(p, a, b, N)
        if val == 0
            # no root
            continue
        elseif val == 1
            J.N[], J.Depth[] = N, depth′             
            push!(Internal, J)
            continue
        end


        setprecision(L) do
            w = b-a
            m = a + (b-a)/2
            c = admissiblepoint(p, m, w/32, n′)
            N′ = max(4, ceil(Int, sqrt(N)))
            push!(Internal, 𝑰(a,c, N′, depth′))
            push!(Internal, 𝑰(c,b, N′, depth′))
        end
    end
    
    return State(Isol, Unresolved, NTuple{n+1,eltype(p)}(p))
    

end


