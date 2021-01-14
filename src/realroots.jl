
## Find precision of a big float
precision(x::Float64) = 53
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


## The Taylor shift here is the most expensive operation
## https://arxiv.org/pdf/1605.00410.pdf has a better strategy
## of partial Taylor shifts, using just the nearby roots
##
## Some comments [here](https://math.stackexchange.com/questions/694565/polynomial-shift)
## are interesting
#
# compute [p(a), p'(a), 1/2p''(a). 1/3! p'''(a)...]
# as q(x) = p(x+a) = p(a) + p'(a)x + 1/2! p''(a)x^2 + 1/3! p'''(a) x^3 ....
# this uses O(n^2) Horner scheme
" `p(x + λ)`: Translate polynomial left by λ "
function poly_translate!(p, λ::S=1) where {S}
    T = promote_type(eltype(p[1]), S)
    n = length(p)
    dps = zeros(T, n)
    for i in n:-1:1
        for j in n:-1:2
            dps[j] = muladd(dps[j], λ, dps[j-1])
        end
        dps[1] = muladd(dps[1], λ, p[i])
    end

    p .= dps
    nothing

end
Tλ(p, lambda=1)   = poly_translate!(copy(p), lambda)



" `R(p)` finds  `x^n p(1/x)` which is a reversal of coefficients "
function poly_reverse!(p)
    reverse!(p)
end
R(p) = poly_reverse!(copy(p))

" `p(λ x)`: scale x axis by λ "
function poly_scale!(p, lambda)
    for i in 2:length(p)
        p[i] *= lambda^(i-1)
    end
end
Hλ(p::Vector{T}, λ=one(T)) where {T} = poly_scale!(copy(p), λ)

"q = p(-x)"
function poly_flip!(p)
    for i in 2:2:length(p)
        p[i] = -p[i]
    end
end


# shift polynomial so (a,b) -> (0, oo)
function poly_shift(p::Vector{T}, a, b) where {T}
    p1 = copy(p)
    poly_translate!(p1, a)
    poly_scale!(p1, b-a)
    poly_reverse!(p1)
    poly_translate!(p1, one(T))
    p1
end


function descartesbound(p, a, b, L =  maximum(precision, (a,b)))
    T = BigFloat
    setprecision(L) do
        q = poly_shift(collect(T, p), T(a), T(b))
        descartescount(q, 1/T(2)^L)
    end
end
descartesbound(p) = descartescount(p, 0.0)

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

function admissiblepoint(p, m, ϵ, n)

    T = BigFloat
    L = precision(m)    
    n′ = sum(divrem(n,2))

    for i in 1:4
        tol = 4/big(2.0)^L
        out = setprecision(L) do
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
    return nothing
end

struct Interval{T}
    a::T
    b::T
    N::Base.Ref{Int}
    bnd::Base.Ref{Int}
end
Interval(a,b, n::Int=4, bnd::Int=-2) = Interval(big(a), big(b), Ref(n), Ref(bnd))

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

struct State
    Internal::Vector{Interval{BigFloat}}                    # DesBound > 1
    Isol::Vector{Interval{BigFloat}}                        # DesBound == 1
    Unresolved::Vector{Interval{BigFloat}}
    p
end

State(p)  = State(BigFloat[], BigFloat[], BigFloat[], ntuple(i -> p[i], length(p)))
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
    println(io, "There were $(length(Is)) isolating intervals found:")
    
    for I in Is
        println(io, I)
    end
end


function tₐ(pa)
    apa = abs(pa)
    L = precision(pa)
    apa > 1/big(2)^L && return ceil(Int, log2(apa))
    return 999999
end

zerotest(st, I) = zerotest(st.p, I.a, I.b)
function zerotest(p, a, b)

    L′ = maximum(precision, (a,b))
    pa, pb = setprecision(L′) do
        evalpoly(a, p), evalpoly(b, p)
    end
    ta, tb = tₐ(a), tₐ(b)

    n = length(p) + 1
    n′ = sum(divrem(n,2))

    L = max(24, L′, max(1, -min(ta-1, tb-1) + 2(n+1) + 1))

    mI = a + (b-a)/2

    bnda = descartesbound(p, a, mI, L)
    bndb = descartesbound(p, mI, b, L)

    iszero(bnda) && iszero(bndb) && return true
    return false
end


onetest(p, a, b) = onetest(State(p), Interval(a,b))

function onetest(st, I)

    p = st.p
    a, b= I

    L′ = maximum(precision, (a,b))
    pa, pb = setprecision(L′) do
        evalpoly(a, p), evalpoly(b, p)
    end
    ta, tb = tₐ(a), tₐ(b)

    n = length(p) + 1
    n′ = sum(divrem(n,2))

    ϵ = (b-a) / (4n)

    ## XXX admissible
    mstar = admissiblepoint(p, a + (b-a)/2, ϵ, n′)
    pmstar = evalpoly(mstar, p)
    t = tₐ(mstar)
    
    L = max(24, L′, max(1, -min(ta-1, tb-1, t-1) + 4(n+2)))
        
    bnda = descartesbound(p, a, mstar)
    bndb = descartesbound(p, mstar, b)
    
    iszero(bnda) && isone(bndb) && return true
    isone(bnda) && iszero(bndb) && return true

    ## update bnd...
    
    return false

end

function newtontest(st, I)

    p = st.p
    a, b = I
    m, w  = a + (b-a)/2, b-a
    
    N::Int = I.N[]
    bnd::Int = I.bnd[]

    p, p′ = p, poly_deriv(p) 
    n = length(p) + 1
    ϵ = 1/(2^ceil(Int, 5 + log(n))) 

    wI = (b-a)
    ξs = [admissiblepoint(p, a + i*wI/4, ϵ * wI, n) for i in 1:3]
    for ξ in ξs
        ξ === nothing && return (false, I)
    end
    
    vs = (evalpoly(ξs[1],p)/evalpoly(ξs[1],p′),
          evalpoly(ξs[2],p)/evalpoly(ξs[2],p′),
          evalpoly(ξs[3],p)/evalpoly(ξs[3],p′))
    L = 53
    ctr = 0

    for _ in 1:10
        ## XXX unclear where L goes here...
        for i in (1,2,3)
            for j in (i+1):3
                # (25)
                abs(vs[i]) > wI && abs(vs[j]) > wI && continue
                abs(vs[i] - vs[j]) < wI/(4n) && continue
                
                if abs(vs[i]) < 2*wI && abs(vs[j]) < 2*wI && abs(vs[i]-vs[j]) > wI/(8n)
                    k̄ = (ξs[j] - ξs[i]) / (vs[j] - vs[i])
                    λᵢⱼ = ξs[i] - k̄ * vs[i] # k cluster is supposed
                    a <= λᵢⱼ <= b  || continue
                    lᵢⱼ = floor(Int, Float64((λᵢⱼ - a) * (4 * N) / wI))

                    aᵢⱼ = a + max(0,  lᵢⱼ - 1) * wI / (4N)
                    bᵢⱼ = a + min(4N, (lᵢⱼ + 2)) * wI/(4N)

                    a⁺ᵢⱼ = (aᵢⱼ == a) ? aᵢⱼ : admissiblepoint(p, aᵢⱼ, ϵ * wI/(4N), n)
                    (!zerotest(st, Interval(a, a⁺ᵢⱼ))) && continue
                    b⁺ᵢⱼ = (bᵢⱼ == a) ? bᵢⱼ : admissiblepoint(p, bᵢⱼ, ϵ * wI/(4N), n)
                    (!zerotest(st, Interval(b⁺ᵢⱼ, b)) ) && continue                    
                    return true, Interval(a⁺ᵢⱼ, b⁺ᵢⱼ, N*2, bnd)

                end
            end
        end
    end

    return (false, I)

end

function boundarytest(st, I)
    T = BigFloat
    a, b = I
    w = b-a
    m = a + w/2

    p = st.p
    n = length(p)-1
    n′ = sum(divrem(n,2))

    N::Int = I.N[]

    wI = b - a
    mₗ, mᵣ = a + wI/(2N), b - wI/(2N)
    ϵ = 1/T(2)^(2 + ceil(Int, log(n)))

    mₗ⁺, mᵣ⁺ = admissiblepoint(p, mₗ, ϵ*wI/N, n′), admissiblepoint(p, mᵣ, ϵ*wI/N, n′)
    a < mₗ⁺ <= mᵣ⁺ < b || return (false, I)

    Iₗ, Iᵣ = Interval(mₗ⁺, b), Interval(a, mᵣ⁺)
    zₗ, zᵣ = zerotest(st, Iₗ), zerotest(st, Iᵣ)

    zₗ && zᵣ && return (true, Interval(mₗ⁺, mᵣ⁺, N, 0))
    zₗ && return (true, Interval(a, mₗ⁺))
    zᵣ && return (true, Interval(mᵣ⁺, b))
    (false, I)

end

        
 ## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function upperbound(p)
    T = BigFloat
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

function lowerbound(p)
    q = p[findfirst(!iszero, p):end]
    poly_flip!(q)
    -upperbound(q)
end
   
"""
    ANewDsc(p; [m=lowerbound(p)], [M=upperbound(p)])
    real_roots_sqfree(p; [m=], [M=])
    real_roots(p; [m=], [M=])

A method to find isolating intervals for the real roots of the polynomial specified by `p`.

* `p`: the polynomial coefficients, `[a₀, a₁, …, aₙ]`, of a **square-free** polynomial.
* `m`: a lower bound for the smallest possible real root
* `M`: an upper bound for the largest possible real root

Returns a `State` instance which has components:

* `Isol` holding the isolating intervals. Iteration over a `State` object will iterate over `Isol`.
* `Unresolved` holding any unresolved intervals. The show method alerts the presence of any such intervals.


Examples:

```
julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1] # mignotte polynomial with two nearby roots

julia> st = ANewDsc(ps)
There were 3 isolating intervals found:
[6.0…, 10.5…]₂₅₆
[0.00787401576326…, 0.00787401637058…]₂₅₆
[0.00787401515549…, 0.00787401576326…]₂₅₆

julia> ps = [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]; # πᵢ₌₁¹⁰ (x-i)

julia> ANewDsc(ps)
There were 10 isolating intervals found:
[9.75…, 10.5…]₂₅₆
[8.5…, 9.5…]₂₅₆
[7.5…, 8.5…]₂₅₆
[6.5…, 7.5…]₂₅₆
[5.75…, 6.75…]₂₅₆
[4.5…, 5.75…]₂₅₆
[3.38…, 4.5…]₂₅₆
[2.5…, 3.38…]₂₅₆
[1.31…, 2.5…]₂₅₆
[0.625…, 1.31…]₂₅₆

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

julia> ANewDsc(ps, m=-4.0, M=0.0)
There were 15 isolating intervals found:
[-0.028107…, 0.0…]₂₅₆
[-0.1123…, -0.02808…]₂₅₆
[-0.2578…, -0.1123…]₂₅₆
[...]
```

To find the zeros, the `Roots` package may be used. For example:

```
julia> using Roots

julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1];

julia> st = ANewDsc(ps)
There were 3 isolating intervals found:
[6.0…, 10.5…]₂₅₆
[0.00787401576326…, 0.00787401637058…]₂₅₆
[0.00787401515549…, 0.00787401576326…]₂₅₆

julia>  [find_zero(st, (I.a, I.b), Roots.BisectionExact()) for I in st]
3-element Array{BigFloat,1}:
 6.93943740962139212443671349244761027220068050171218581650766763204507611476287
 0.007874016089132754403608727898779727134193464194254785228308443233139694038151204
 0.007874015406930341157555003028161633376551552518768059431667490175426147684590201
```

The default bracketing method for `BigFloat` can have issues with some problems, so we use `BisectionExact` above. Alternatively, `Roots.bisection(st, I...)` may be used.

The algorithm used is a simplification of

Computing Real Roots of Real Polynomials ... and now For Real!
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
arXiv:1605.00410; DOI:	10.1145/2930889.2930937

!!! note
    a square free polynomial can be found through `p/gcd(p, p')`, though in practice this
    calculation is numerically unstable.

!!! note
    This is much slower as the function provided through `arblib` in `Hecke.jl`,
    which itself says is not competitive with specialized algorithms, such as provided in the
    RS library of the paper authors.

"""
function ANewDsc(p; m=lowerbound(p), M=upperbound(p))
    DEBUG = false

    st = State(p)
    n = length(p) - 1
    n′ = sum(divrem(n,2))

    I = Interval(promote(float(m), float(M))...)
    push!(st.Internal, I)

    while !isempty(st.Internal)
        I = pop!(st.Internal)

        if zerotest(st, I)
            continue
        end

        if onetest(st, I)
            push!(st.Isol, I)
            continue
        end


        val, J = newtontest(st, I)
        if val
            DEBUG && @show :newton
            push!(st.Internal, J)
            continue
        end
        
        val, J = boundarytest(st, I)
        if val
            DEBUG && @show :bound, I.a, I.b
            push!(st.Internal, J)
            continue
        end

        a,b = I
        pa, pb = setprecision(maximum(precision, (a,b))) do
            evalpoly(a, p), evalpoly(b, p)
        end

        if iszero(pa) || iszero(pb)
            DEBUG && @show :zero_on_boundary
            push!(st.Unresolved, I)
            continue
        end
        
        w = b-a
        m = a + (b-a)/2
        DEBUG && @show :linear
        c = admissiblepoint(p, m, w/4, n′)
        push!(st.Internal, Interval(a,c))
        push!(st.Internal, Interval(c,b))
    end
    
    st
end

## polynomials
# real_roots_sqfree(p; kwargs...) = ANewDsc(p.coeffs; kwargs...)

# function real_roots(p; kwargs...)
#     u,v,w,_,_ = ngcd(p, derivative(p))
#     real_roots_sqfree(v; kwargs...)
# end

# export real_roots
