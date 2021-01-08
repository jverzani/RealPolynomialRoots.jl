## Find real roots of a polynomial using a method from
## "Computing Real Roots of Real Polynomials ... and now For Real!"
## by Alexander Kobel, Fabrice Rouillier, Michael Rouillier
## (https://arxiv.org/abs/1605.00410)
## DOI:	10.1145/2930889.2930937
## 
## Earlier work of a Descartes-like method  [Vincent's theorem](http://en.wikipedia.org/wiki/Vincent's_theorem)
## are here "Efficient isolation of polynomial’s real roots" by
## Fabrice Rouillier; Paul Zimmermann
##
##
## The implementation in Hecke.jl uses the full power of arblib and is
## much better engineered than this.
## This uses the ArbNumerics.jl package
##




## We consider intervals [a,b] with a,b of type ArbFloat{L}
## Interval with [a,b], N
## N is for Newton-Test
struct Interval{T}
    a::T
    b::T
    N::Base.Ref{Int}
    bnd::Base.Ref{Int}
end
Interval(a,b, n::Int=4, bnd::Int=-2) = Interval(_promote(a,b)...,Ref(n),Ref(bnd))
width(I::Interval) = I.b-I.a
midpoint(I::Interval) = I.a + (0.5) * width(I)


## how we store the state of our algorithm to find zero in [0,1]
# we use a tuple to store `p`
# the Intervals may have different precision so this is all untyped
struct State
    Internal                    # DesBound > 1
    Isol                        # DesBound == 1
    Unresolved
    p
end

State(p)  = State(Any[], Any[], Any[], ntuple(i -> p[i], length(p)))
(st::State)(x) = evalpoly(x, st.p)

# **Much** faster to define poly operations over vectors than to
# use `Polynomials.ImmutablePolynomial`.

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

" `p(x + λ)`: Translate polynomial left by λ "
function poly_translate!(p::Vector{T}, lambda=1) where {T}
    p1 = copy(p)
    p[1] = evalpoly(lambda, p1)
    m = one(T)
    for k in 2:length(p1)
        p[k] = evalpoly(lambda, poly_deriv!(p1)) / m
        m *= k
    end
    p
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

## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function _upperbound(p)
    T = ArbFloat{64}    
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

    T = ArbFloat{64}
    _promote(one(T), out)[2]
end

function _lowerbound(p::Vector{T}) where {T <: Real}
    q = p[findfirst(!iszero, p):end]
    
    poly_flip!(q)
    ret = -_upperbound(q)
    ret
end




## Descartes Bounds
sgn(p) = sign(p)

function sgn(p::T) where {L, T <: ArbFloat{L}}
    τ = 1/(T(2)^L)
    p <= -τ && return -1
    p >=  τ && return 1
    return nothing
end



function descartesbound(st::State, I::Interval{T}) where {T}
    bnd = I.bnd[]
    bnd != -2 && return bnd

    bnd = descartesbound(collect(T, st.p), T(I.a), T(I.b))
    I.bnd[] = bnd
    bnd
end

function descartesbound(p, a, b)
    q = poly_shift(p, a, b)
    descartesbound(q)
end

function descartesbound(q)
    cnt = -1
    s = 0
    for qᵢ in q
        sᵢ = sgn(qᵢ)
        sᵢ === nothing && return -1
        if sᵢ != s
            s = sᵢ
            cnt += 1
        end
    end
    return cnt
end
    

# identify 2^(tₐ-1) <= y < t^(tₐ)
tₐ(y) = ceil(Int, log2(abs(y)))

# m[ϵ] evenly spaced points in B(x,ϵ)
function multipoint(m, ϵ, n′)
    (m + i * ϵ/n′ for i in -n′:n′)
end
    

# within B(m, ϵ) find a point which is bigger than a fraction
# of the maximum over the multipoing with a precision of a certain size.
function admissiblepoint(p, m::ArbFloat{LL}, ϵ, n) where {LL}

    n′ = ceil(Int, n/2)
    L = LL
    for i in 1:10
        T = ArbFloat{L}
        tol = 4/T(2.0)^L 
        mx, vx = zero(T), zero(T)
        for mᵢ in multipoint(T(m), T(ϵ), n)
            vᵢ = abs(evalpoly(mᵢ, p))
            if vᵢ > vx
                vx, mx = vᵢ, mᵢ
            end
        end
        vx > tol && return T(mx)
        L = 2 * L
    end
    return nothing
end

function zerotest(st::State, I::Interval{T₀}, R=T₀) where {T₀}

    iszero(I.bnd[]) && !isbracket(st, I) && return true
    
    a,b = I.a, I.b
    p = R == T₀ ? st.p : R(st.p)
    pa, pb = evalpoly(R(a), p), evalpoly(R(b), p)
    ta,tb = tₐ(pa), tₐ(pb)

    n = length(p) + 1
    n′ = sum(divrem(n,2))
    
    L = max(24, arbL(a), max(1, -min(ta-1, tb-1) + 2(n+1)+1))
    ## really need L' here? (Corollary 15)
    wI = b - a
    mI = a + wI/2
    mIa = admissiblepoint(p, mI, wI/8, n′)
    T = ArbFloat{2^(ceil(Int, log2(L)))}

    q = collect(T,p)
    dl = descartesbound(q, T(a), T(mI))
    dr = descartesbound(q, T(mI), T(b))
    if dl >= 0 && dr >= 0
        I.bnd[] = dl + dr
    end
    
    return iszero(dl) && iszero(dr) 
end

isbracket(st::State, I) = isbracket(st.p, I)
isbracket(p, I) = sign(evalpoly(I.a,p)) * sign(evalpoly(I.b,p)) <= 0

function onetest(st, I::Interval{T₀}, R=T₀)  where {T₀}

    p = R == T₀ ? st.p : R.(st.p)    
    a, b = R(I.a), R.(I.b)
    ta,tb = tₐ(evalpoly(a,p)), tₐ(evalpoly(b, p))
    wI = b - a
    n = length(p) + 1
    ϵ = wI/(4n)

    
    mstar = admissiblepoint(p, a + (b-a)/2, ϵ, n)
    t = tₐ(evalpoly(mstar, p))

    L = max(24,  arbL(a), max(1, -min(ta-1, tb-1, t-1) + 4n+2)) 
    T = ArbFloat{2^(ceil(Int, log2(L)))}

    Ia, Ib = Interval(T(a), T(mstar)), Interval(T(mstar), T(b))
    zta, ztb = zerotest(st, Ia, T), zerotest(st, Ib, T)
    # XXX isbracket needed here
    zta && Ib.bnd[] == 1 && isbracket(p, Ib) && return (true, Ib)
    ztb && Ia.bnd[] == 1 && isbracket(p, Ia) && return (true, Ia)    
#    zta && Ib.bnd[] == 1 && return (true, Ib)
#    ztb && Ia.bnd[] == 1 && return (true, Ia)    
    
    return (false, Interval(a, b, I.N[], Ia.bnd[] + Ib.bnd[]))
end

function boundarytest(st, node::Interval{T}) where {T}
    a, b, m, w  = node.a, node.b, midpoint(node), width(node)
    p = st.p
    n = length(p)-1
    n′ = sum(divrem(n,2))
    N::Int = node.N[]
    wI = b - a
    mₗ, mᵣ = a + wI/(2N), b - wI/(2N)
    ϵ = 1/T(2)^(2 + ceil(Int, log(n)))

    mₗ⁺, mᵣ⁺ = admissiblepoint(p, mₗ, ϵ*wI/N, n′), admissiblepoint(p, mᵣ, ϵ*wI/N, n′)
    a < mₗ⁺ <= mᵣ⁺ < b || return (false, node)

    Iₗ, Iᵣ = Interval(mₗ⁺, b), Interval(a, mᵣ⁺)
    zₗ, zᵣ = zerotest(st, Iₗ), zerotest(st, Iᵣ)

    zₗ && zᵣ && return (true, Interval(mₗ⁺, mᵣ⁺, N, 0))
    zₗ && return (true, Interval(a, mₗ⁺))
    zᵣ && return (true, Interval(mᵣ⁺, b))
    (false, node)

end


function newtontest(st, node)

    a, b, m, w  = node.a, node.b, midpoint(node), width(node)
    N::Int = node.N[]
    bnd::Int = node.bnd[]
    p, p′ = st.p, poly_deriv(st.p) 
    n = length(p) + 1
    ϵ = 1/(2^ceil(Int, 5 + log(n))) 

    wI = (b-a)
    ξs = [admissiblepoint(p, a + i*wI/4, ϵ * wI, n) for i in 1:3]
    for ξ in ξs
        ξ === nothing && return (false, node)
    end
    
    vs = (evalpoly(ξs[1],p)/evalpoly(ξs[1],p′),
          evalpoly(ξs[2],p)/evalpoly(ξs[2],p′),
          evalpoly(ξs[3],p)/evalpoly(ξs[3],p′))
    L = 64
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

    return (false, node)

end


function linearstep(st, node::Interval{T}) where {L, T<:ArbFloat{L}}
    a, b, N, bnd = node.a, node.b, node.N[], node.bnd[]
    p = st.p
    n = length(p) - 1
    n′ = sum(divrem(n, 2))
    N′ = max(4, floor(Int, sqrt(N)))
    ϵ = (b-a)/T(2)^(ceil(Int, 2 + log(n)))
    mstar = admissiblepoint(collect(T,p), a + (b-a)/2, ϵ, n′)
    I1, I2 = Interval(a, mstar, N′), Interval(mstar, b, N′)    
    I1, I2
end


function nextstep!(st::State, node) 

    if zerotest(st, node)
        if isbracket(st, node)
#            @show :this_shouldnot_happen, node.bnd[]
            push!(st.Unresolved, node)
        end
        return nothing
    end
    
    val, J = onetest(st, node)
    if val
        push!(st.Isol, J)
        return nothing
    end
    
    for test in (newtontest,  boundarytest)
        val, I = test(st, node)
        if val
#            @show test
            push!(st.Internal, I)
            return nothing
        end
    end

#    @show :linear
    for I in linearstep(st, node)
        bnd = descartesbound(st, I)
        if bnd == -1
            push!(st.Unresolved, I)
        else
            push!(st.Internal, I)
        end
    end
    return nothing
end

# assume p is square free
function anewdsc(p; m=_lowerbound(p), M=_upperbound(p))
    st = State(p)
    I = Interval(m, M)
    push!(st.Internal, I)
    while !isempty(st.Internal)
        I = pop!(st.Internal)
        nextstep!(st, I)
    end
    st
    
end

"""
    real_roots_sqrfree(p; [m=], [M=])

Call `anewdsc` to find isolating intervals, then find zeros using a bracketing method

`p`: the polynomial coefficients, `[a₀, a₁, …, aₙ]`, of **square-free** polynomial.
`m`, `M`: used to narrow search space 

Examples:
```
julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1] # mignotte polynomial with two nearby roots

julia> real_roots_sqfree(ps)
3-element Array{ArbFloat,1}:
 6.93943740962
 0.00787401608913275440360872789878
 0.007874015406930341157555003028162

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

julia> real_roots_sqfree(ps, m=-4, M=0)
15-element Array{ArbFloat{64},1}:
 -0.00701539819528
 -0.0629983810833
 -0.174209253846
 -0.339102025104
 -0.555265632219
 -0.819305999802
 -1.12667548938
 -1.47143642655
 -1.84593922753
 -2.24038666835
 -2.64224182084
 -3.03541295495
 -3.39910134596
 -3.70611151706
 -3.92061671115
```

The algorithm used is based on 

Computing Real Roots of Real Polynomials ... and now For Real!
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
arXiv:1605.00410; DOI:	10.1145/2930889.2930937

!!! note
    a square free polynomial can be found through `p/gcd(p, p')`, though in practice this
    calculation is numerically unstable.

!!! note
    This is not nearly as fast as the function provided through `arblib` in `Hecke.jl`.

"""
function real_roots_sqfree(p; kwargs...)

    st = anewdsc(p; kwargs...)

    out = [find_zero(st, (I.a, I.b)) for I in st.Isol]
    if !isempty(st.Unresolved)
        println("Some intervals found were unresolved")
        for I in st.Unresolved
            println(I)
        end
    end

    out
end
    

## --------------------------------------------------

# """

#      real_roots(p, [m], [M]; square_free=true)

# Returns real roots of a square-free polynomial presented via its coefficients
# `[p_0, p_1, ..., p_n]`. 

# * `p`: polynomial coefficients, `Vector{T<:Real}`
# * `m`: lower bound on real roots. Defaults to `_lowerbound(p)`
# * `M`: upper bound on real roots. Defaults to `_upperbound(p)`
# * `square_free`::Bool. If false, the polynomial `agcd(p, polyder(p))` is used. This polynomial---in theory--- would have the
# same real roots as `p`, however in practice the approximate `gcd` can be off.

# """
# real_roots(p::Poly{T}, args...; kwargs...) where {T <: Real} = real_roots(p.a, args...; kwargs...)
# function real_roots(p::Vector{T}, m = lowerbound(p), M=upperbound(p); square_free::Bool=true) where {T <: Real}

#     # deflate zero
    
#     nzroots = findfirst(!iszero, p) - 1
#     if nzroots > 0
#         p = p[1+nzroots:end]
#     end

#     if !square_free
#         error("not implemented ")
#     end

    
#     st = isolate_roots(p, m, M)
    
#     if length(st.Unresolved) > 0
#         println("Some intervals are unresolved:")
#         println("------------------------------")
#         for node in st.Unresolved
#             @show node
#             #@printf "* There may be up to %d roots in (%0.16f, %0.16f).\n" node.bnd[] node.a node.b
#         end
#         println("------------------------------")                
#     end

#          rts = zeros(T, length(st.Isol))
#      for i in eachindex(st.Isol)
#          node = st.Isol[i]
#          a, b = node.a, node.b

         
#          pa, pb = evalpoly(a, p), evalpoly(b, p)
#          if sign(pa) * sign(pb) >= 0
#              rt = abs(pa) > abs(pb) ? b : a
#          else
#              # for higher precision than Float64, these no longer are guaranteed to converge...
#              rt = NaN*one(T)
#              for P in (Roots.A42(), Roots.AlefeldPotraShi(), Roots.BisectionExact())
                 
#                  zp = Roots.ZeroProblem(P, x->evalpoly(x,p), (a,b))                 
#                  for _ in zp; end
#                  rt = Roots.decide_convergence(zp)
#                  !isnan(rt) && break
#              end
#          end

#          rts[i] = rt
#      end

#     nzroots > 0 && push!(rts, zero(T))
#     rts
# end
        
