

## Find real roots of a polynomial using a DesCartes Method
##
## State of the art for this approach is summarized here:
## Cf. "Computing Real Roots of Real Polynomials ... and now For Real!" by Alexander Kobel, Fabrice Rouillier, Michael Rouillier
## (https://arxiv.org/abs/1605.00410) 
## 
## Earlier work of a Descartes-like method  [Vincent's theorem](http://en.wikipedia.org/wiki/Vincent's_theorem)
## are here "Efficient isolation of polynomial’s real roots" by
## Fabrice Rouillier; Paul Zimmermann
##
## This implementation doesn't take nearly enough care with the details, but takes some ideas
## to implement a means to find real roots of non-pathological polynomials (lowish degree, roots separated)
##
## XXX Needs work XXX

## Polynomial transformations
##
## The Taylor shift here is the most expensive operation
## https://arxiv.org/pdf/1605.00410.pdf has a better strategy
## of partial Taylor shifts, using just the nearby roots
## 


#Polynomials.degree(p::Vector{T}) where {T} = findlast(!iszero,p) - 1

_promote(x,y) = promote(x,y)
function _promote(x::Float64, y::Float64)
    (ArbFloat{64}(x),ArbFloat{64}(y))
end
function _promote(x::ArbFloat{L}, y::Float64) where {L}
    (x, ArbFloat{L}(y))
end

function _promote(x::Float64, y::ArbFloat{L}) where {L}
    return (ArbFloat{L}(x),y)
end

function _promote(x::ArbFloat{L}, y::ArbFloat{M}) where {L, M}
    L > M && return (x, ArbFloat{L}(y))
    return (ArbFloat{M}(x), y)
end

function _promote(x::ArbFloat, y::BigFloat)
    (big(x), y)
end
    
arbL(x::ArbFloat{L}) where {L} = L


## Interval with [a,b], N
## N is for Newton-Test
struct Interval{T,S}
    a::T
    b::S
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

# poly is [p_0, p_1, ..., p_n] where there may be zeros in last terms

# p + q
function poly_add!(p::Vector{T}, q::Vector{S}) where {T, S}
    n,m = length(p), length(q)
    l = min(n, m)
    for i in 1:l
        p[i] += q[i]
    end
    for i in l:m
        push!(p, q[i])
    end
end

# p-q
poly_sub!(p::Vector{T}, q::Vector{S}) where {T, S} = poly_add(p, -q)

# may be more than one way
function poly_mul!(p1::Vector{T}, p2::Vector{S}) where {T, S}
    R = promote_type(T,S)

    n = length(p1)
    m = length(p2)
    
    a = zeros(R,m+n+1)

    for i = 1:n
        for j = 1:m
            a[i+j+1] += p1[i] * p2[j]
        end
    end
    p1[:] = a

end


function poly_deriv!(p::Vector{T}) where {T}
    for i = 1:length(p)-1
        p[i] = i * p[i+1]
    end
    p[end] = zero(T)
    p
end



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
    #    p0 = Tλ(p, a)
    #    p1 = Hλ(p0, b-a)
    #    p2 = Tλ(R(p1),1)

    p1 = copy(p)
    poly_translate!(p1, a)
    poly_scale!(p1, b-a)
    poly_reverse!(p1)
    poly_translate!(p1, one(T))
    p1
end

## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function upperbound(p)
    T = ArbFloat{64}
    p = p[findfirst(!iszero, p):end]
    descartesbound(collect(T,p)) == 0 && return zero(T)

    p = p[findfirst(!iszero, p):end]
    
    q, d = p/p[end], length(p)-1
    
    d == 0 && error("degree 0 is a constant")
    d == 1 && abs(q[1])


    a1 = abs(q[d])
    B = maximum([abs(q[i]) for i in 1:(d-1)])

    a,b,c = 1, -(1+a1), a1-B
    out = (-b + sqrt(b^2 - 4a*c))/2
    T(out)
end

function lowerbound(p::Vector{T}) where {T <: Real}
    q = p[findfirst(!iszero, p):end]
    
    poly_flip!(q)
    ret = -upperbound(p)
    ret
end





## Descartes Bounds

" Descartes bound on (0,oo). Just count sign changes"
# function descartes_bound(p::Vector{T}) where {T}
#     length(p) == 0 && return -1
#     cnt, sgn = 0, sign(p[1])
#     for i in 2:length(p)
#         nsgn = sign(p[i])
#         if nsgn * sgn < 0
#             sgn = nsgn
#             cnt += 1
#         end
#     end
#     cnt
# end

# function sgn(I::IntervalArithmetic.Interval{T}) where {T}
#     z = zero(T)
#     z < I && return 1
#     z > I && return -1
#     (I.hi-I.lo) <= 1000*sqrt(eps(T)) && return 0
#     return nothing #  sentinel
# end
sgn(p) = sign(p)

function sgn(p::T) where {L, T <: ArbFloat{L}}
    τ = 100/(T(2)^L)
    p <= -τ && return -1
    p >=  τ && return 1
    return nothing
end


function descartes_count(ps)
    cnt = -1
    s = 0
    for pᵢ in ps
        sᵢ = sgn(pᵢ)
        #        @show pᵢ, sᵢ, eps(eltype(pᵢ))
        sᵢ === nothing && error("increase precision")
        iszero(sᵢ) && continue
        if sᵢ != s
            s = sᵢ
            cnt += 1
        end
    end
    cnt
end

function descartesbound(st::State, I::Interval{T,T}) where {T}
    bnd = I.bnd[]
    bnd != -2 && return bnd

    
#    bnd = descartesbound(collect(T, st.p), T(I.a), T(I.b))  ## DOUBLE XXX
    
    L = arbL(I.a)
    TT = ArbFloat{2L}
    
    bnd = descartesbound(collect(TT, st.p), TT(I.a), TT(I.b))  ## DOUBLE XXX
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
    

# split interval if not precise enough?
# function descartes_bound(p::Vector{T}, S=T) where {T}
#     ϵ = sqrt(eps(S))
#     qs = [IntervalArithmetic.Interval(pᵢ-ϵ, pᵢ+ϵ) for pᵢ in p]
#     descartes_count(qs)
# end
# function descartes_bound(p::Vector{T}, a, b, S=T) where {T}

#     # should not have degree(p) <= 1
#     ϵ = sqrt(eps(S))
#     qs = [IntervalArithmetic.Interval(pᵢ-ϵ, pᵢ+ϵ) for pᵢ in p]
#     descartes_count(poly_shift(qs, a, b))
# end

" Descartes bound on (a, b)"
#descartes_bound(p::Vector{T}, a, b) where {T} = descartes_bound(poly_shift(p, a, b))
function DescartesBound(st, node)
    bnd::Int = node.bnd[]
    bnd != -2 && return bnd # already done if not -2
    
    bnd = descartes_bound(st.p, node.a, node.b)
    node.bnd[] = bnd
    bnd
end

## Tests

## Return true or false
zero_test(st::State, node)::Bool  =   DescartesBound(st, node) == 0 
one_test(st::State, node)::Bool  =   DescartesBound(st, node) == 1  

## return count -1 (can't tell), 0, 1, or more
zero_one_test(st::State, node) = DescartesBound(st, node)  


# find admissible point
# XXX improve me
# function find_admissible_point(st::State{T},  I::Interval, m=midpoint(I), Ni::T=one(T), c::T=one(T)) where {T}
#     N = ceil(Int, Float64(Polynomials.degree(st.p)/2))
#     ep = min(m-I.a, I.b - m) / (4*Ni)
#     mis = [m + i/N * ep for i in -N:N]
#     curmin = min(norm(st(I.a)), norm(st(I.b)))/100
#     for m in mis #shuffle(mis)
#         (m < I.b || m > I.a) || continue
#         descartes_bound(st.p, I.a, m) == -1 && continue
#         descartes_bound(st.p, m, I.b) == -1 && continue        
#         norm(st(m)) > 0 && return m
#     end
#     @show I
#     error("No admissible point found")
# #    mx, i = findmax(norm.(st.p.(mis)))
# #    mis[i]
# end

# n = degree of P
# n′ = ceil(n/2)
# var(P, I1) + var(P, I2) ≤ var(P, I)

# identify 2^(tₐ-1) <= y < t^(tₐ)
tₐ(y) = ceil(Int, log2(abs(y)))

M(z) = max(1, abs(z))

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
        tol = 2*2*4/T(2.0)^L
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

function zerotest(st::State, I::Interval{T₀,T₀}, R=T₀) where {T₀}

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
    T = ArbFloat{2*2^(ceil(Int, log2(L)))}

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

function onetest(st, I::Interval{T₀,T₀}, R=T₀)  where {T₀}

    p = R == T₀ ? st.p : R.(st.p)    
    a, b = R(I.a), R.(I.b)
    ta,tb = tₐ(evalpoly(a,p)), tₐ(evalpoly(b, p))
    wI = b - a
    n = length(p) + 1
    ϵ = wI/(4n)

    
    mstar = admissiblepoint(p, a + (b-a)/2, ϵ, n)
    t = tₐ(evalpoly(mstar, p))

    L = 2max(24,  arbL(a), max(1, -min(ta-1, tb-1, t-1) + 4n+2)) ## need 2 here for mignotte(31,14)
    T = ArbFloat{2^(ceil(Int, log2(L)))}

    Ia, Ib = Interval(T(a), T(mstar)), Interval(T(mstar), T(b))
    zta, ztb = zerotest(st, Ia, T) && !isbracket(st, Ia), zerotest(st, Ib, T) && !isbracket(st, Ib)
    zta && Ib.bnd[] == 1 && isbracket(p, Ib) && return (true, Ib)
    ztb && Ia.bnd[] == 1 && isbracket(p, Ia) && return (true, Ia)    
    
    return (false, Interval(a, b, I.N[], Ia.bnd[] + Ib.bnd[]))
end

function boundarytest(st, node)
    a, b, m, w  = node.a, node.b, midpoint(node), width(node)
    p = st.p
    n = length(p)-1
    n′ = sum(divrem(n,2))
    N::Int = node.N[]
    wI = b - a
    mₗ, mᵣ = a + wI/(2N), b - wI/(2N)
    ϵ = 1/2^(2 + ceil(Int, log(n)))

    mₗ⁺, mᵣ⁺ = admissiblepoint(p, mₗ, ϵ*wI/N, n′), admissiblepoint(p, mᵣ, ϵ*wI/N, n′)
    a < mₗ⁺ <= mᵣ⁺ < b || return (false, node)

    Iₗ, Iᵣ = Interval(a, mₗ⁺), Interval(mᵣ⁺, b)
    zₗ, zᵣ = zerotest(st, Iₗ), zerotest(st, Iᵣ)

    ##
    zₗ && isbracket(st, Iₗ) && @show :huh, Iₗ
    zᵣ && isbracket(st, Iᵣ) && @show :huh, Iᵣ
    
    zₗ && zᵣ && return (true, Interval(mₗ⁺, mᵣ⁺, N, node.bnd[]))
    zₗ && return (true, Interval(mₗ⁺, b, N, node.bnd[]))
    zᵣ && return (true, Interval(a, mᵣ⁺, N, node.bnd[]))
    (false, node)

end


function newtontest(st, node)

    a, b, m, w  = node.a, node.b, midpoint(node), width(node)
    N::Int = node.N[]
    bnd::Int = node.bnd[]
    p, p′ = st.p, poly_deriv!(collect(st.p))
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
                #@show vs[i], vs[j], wI
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
                    (!zerotest(st, Interval(a, a⁺ᵢⱼ)) || isbracket(st, Interval(a, a⁺ᵢⱼ))) && continue
                    b⁺ᵢⱼ = (bᵢⱼ == a) ? bᵢⱼ : admissiblepoint(p, bᵢⱼ, ϵ * wI/(4N), n)
                    (!zerotest(st, Interval(b⁺ᵢⱼ, b)) || isbracket(st, Interval(b⁺ᵢⱼ, b))) && continue                    
                    
                    ### XXXreturn true, Interval(a⁺ᵢⱼ, b⁺ᵢⱼ, N^2, bnd)
                    return true, Interval(a⁺ᵢⱼ, b⁺ᵢⱼ, N*2, bnd)

                end
            end
        end
    end

    return (false, node)

end


function linearstep(st, node::Interval{T,T}) where {L, T<:ArbFloat{L}}
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
#-0.853490346223, -0.0131314098584
# find splitting point
# find an admissible point that can be used to split interval. Needs to have all
# coefficients not straddle 0
# return (logical, two intervals)
# function split_interval(st::State{T},I::Interval,  m=midpoint(I), Ni=one(T), c=one(T)) where {T}
#     N = ceil(Int, Float64(Polynomials.degree(st.p)/2))
#     ep = min(1, width(I)) / (16*Ni)
#     mis = T[m + i/N * ep for i in -N:N]
#     mis = filter(m -> m > I.a && m < I.b, mis)
#     mx, i = findmax(norm.(st.(mis)))

#     ## Now, we need a point that is bigger than max and leaves conclusive
#     for i in eachindex(mis)
#         mi = mis[i]
#         abs(st(mi)) >= min(mx/4, one(T)) || continue
#         ileft = Interval(I.a, mi, I.N[], -2)
#         nl = DescartesBound(st, ileft)
#         #nl = descartes_bound(fatten(st.p), ileft.a, ileft.b)                
#         nl == -1 && continue
#         iright = Interval(mi, I.b, I.N[], -2)
#         nr = DescartesBound(st, iright)
#         #nr = descartes_bound(fatten(st.p), iright.a, iright.b)                        
#         nr == -1 && continue

#         # identify degenerate cases here. This is not good.
# #        if nl + nr < I.bnd[]
# #            @warn "possible numeric issue; increase precision?"
# #        end
        
#         ileft.bnd[] = nl; iright.bnd[] = nr

        
#         return (true, ileft, iright)
#     end
#     #    println("DEBUG: $I is a bad interval for splitting?")
#     return (false, I, I)
# end


# ## split interval
# ## adds to intervals if successful
# ## adds node to Unresolved if not
# function linear_step(st::State{T}, node) where {T}
#     succ, I1, I2 = split_interval(st, node)

    
    
#     if succ
#         N::Int = node.N[]
#         N′ = max(4, round(Int, sqrt(N)))
#         for I in (I1, I2)
#             bnd = DescartesBound(st, I)
#             iszero(bnd) && continue
#             if bnd == 1
#                 push!(st.Isol, I)
#             else 
#                 I.N[] = N′
#                 push!(st.Internal, I)
#             end
#         end

# #        b0 = node.bnd[]
# #        b1, b2 = I1.bnd[], I2.bnd[]
# #        if b0 > b1 + b2
# #            @show b0, b1, b2
# #            @show poly_shift(st.p, I1.a, I1.b)
# #            @show poly_shift(st.p, I2.a, I2.b)
# #        end
#     else
#         push!(st.Unresolved, node)
#     end
#     return true
    
# end

function nextstep!(st::State, node) 

    if zerotest(st, node)
        if isbracket(st, node)
            @show :this_shouldnot_happen, node.bnd[]
            push!(st.Unresolved, node)
        end
        return nothing
    end
    
    val, J = onetest(st, node)
    if val
        @show :onetest, node.a, node.b, J.a, J.b
        push!(st.Isol, J)
        return nothing
    end
    
    for test in (newtontest,)#,  boundarytest)
        val, I = test(st, node)
        if val
            @show test
            push!(st.Internal, I)
            return nothing
        end
    end

    @show :linear
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


function anewdsc(p; m=lowerbound(p), M=upperbound(p))
    st = State(p)
    I = Interval(m, M)
    push!(st.Internal, I)

    while !isempty(st.Internal)
        I = pop!(st.Internal)
        #@show :consider, I.a, I.b
        nextstep!(st, I)
    end
    st
end
        

# ## return (true, I), (false, node)
# ## I will be smaller interval containing all roots in node
# function newton_test(st::State{T}, node) where {T}
#     NMAX = 1024  # 2147483648 = 2^31
#     (node.N[] > NMAX) && return (false, node) 

#     a, b, m, w  = node.a, node.b, midpoint(node), width(node)
#     N::Int = node.N[]
#     bnd::Int = node.bnd[]
#     p, p′ = st.p, poly_deriv!(copy(st.p))


#     wI = (b-a)
#     ξs = find_admissible_point.((st,), (node,), a .+ (1,2,3) .* (wI/4))
#     vs = (evalpoly(ξs[1],p)/evalpoly(ξs[1],p′),
#           evalpoly(ξs[2],p)/evalpoly(ξs[2],p′),
#           evalpoly(ξs[3],p)/evalpoly(ξs[3],p′))
#     for i in (1,2,3)
#         for j in (i+1):3
#             k̄ = (ξs[j] - ξs[i]) / (vs[i] - vs[j])
#             λᵢⱼ = ξs[i] + k̄ * vs[i]
#             a <= λᵢⱼ <= b || continue
#             lᵢⱼ = floor(Int, Float64((λᵢⱼ - a) * 4 * N / wI))
#             aᵢⱼ = a + max(0,  lᵢⱼ - 1) * wI / (4N)
#             bᵢⱼ = a + min(1, (lᵢⱼ + 2)/(4N)) * wI
#             a⁺ᵢⱼ = (aᵢⱼ == a) ? aᵢⱼ : find_admissible_point(st, node, aᵢⱼ)
#             b⁺ᵢⱼ = (bᵢⱼ == a) ? bᵢⱼ : find_admissible_point(st, node, bᵢⱼ)
#             Il = Interval(a, a⁺ᵢⱼ , N, -2)
#             Ir = Interval(b⁺ᵢⱼ , b, N, -2)
#             if zero_test(st, Il) &&
#                 zero_test(st, Ir)
#                 return (true, Interval(a⁺ᵢⱼ, b⁺ᵢⱼ, N^2, bnd))
#             end
#         end
#     end

               
#     ## boundary test
#     mlstar::T = find_admissible_point(st, node, a + wI/(2N))
#     if mlstar > a && zero_test(st, Interval(mlstar, b, N, -2))
#         return (true, Interval(a, mlstar, N, bnd))
#     end

#     mrstar::T = find_admissible_point(st, node, b - w/(2N))
#     if mrstar < b && zero_test(st, Interval(a, mrstar, N, -2))
#         return (true, Interval(mrstar, b, N, bnd))
#     end
    
#     return (false, node)
# end

# ## Add successors to I
# ## We have
# function addSucc(st::State{T}, node) where {T}

#     val, I = newton_test(st, node)

#     if val
# #        @show :newton
#         push!(st.Internal, I)
#     else
# #        @show :linear
#         succ = linear_step(st, node)
#         #!succ && @warn("node $node was a failure")
#     end
#     true
# end

# ## m, M should bound the roots
# ## essentially algorithm 4
# function ANewDsc(p::Vector{T}, m = lowerbound(p), M=upperbound(p)) where {T <: Real}

#     st = State(p)
#     base_node = Interval(m, M, 4, -2)    
#     DescartesBound(st, base_node)


#     if base_node.bnd[] == -1
#         append!(st.Internal, break_up_interval(st, base_node, 4))
#     else
#         push!(st.Internal, base_node)
#     end

#     ctr = 0
#     while length(st.Internal) > 0
#         node = pop!(st.Internal)
#         ctr +=1;
# #        @show ctr, node
#         bnd = DescartesBound(st, node) 
#         if bnd < 0
#             push!(st.Unresolved, node)
#             continue
#         elseif bnd == 0
#             continue
#         elseif bnd == 1
#             push!(st.Isol, node)
#             continue
#         else
#             a, b = node.a, node.b
#             if abs(a-b) <= degree(p) * sqrt(eps(T))
#                 push!(st.Unresolved, node)
#                 continue
#             else
#                 addSucc(st, node)
#             end
#         end

#         #@show sum(I.bnd[] for I in st.Internal) + length(st.Isol)   
        
#     end
#     st
# end


# # populate `Isol`
# # p must not have any roots with even degree. (e.g. no (x-c)^(2k) exists as a factor for any c,k
# # assumed square free (at least no roots of even multiplicity)
# function isolate_roots(p::Vector{T}, m, M) where {T <: Real}

# #    try
#     st = ANewDsc(p, m, M)
#     return st
#     # catch err
#     #     if  !(T <: BigFloat)
#     #         try
#     #             st = ANewDsc(convert(Poly{BigFloat}, p), m, M)
#     #             return st
#     #         catch err
#     #             rethrow(err)
#     #         end
#     #     end
#     # end
        
# end



# # for interface
# struct PolyType{T}
#     p::Vector{T}
# end
# (p::PolyType)(x) = evalpoly(x, p.p)



# """

#      real_roots(p, [m], [M]; square_free=true)

# Returns real roots of a polynomial presented via its coefficients `[p_0, p_1, ..., p_n]`. 

# * `p`: polynomial coefficients, `Vector{T<:Real}`
# * `m`: lower bound on real roots. Defaults to `lowerbound(p)`
# * `M`: upper bound on real roots. Defaults to `upperbound(p)`
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
        
