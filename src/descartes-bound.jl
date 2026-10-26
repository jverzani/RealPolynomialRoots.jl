Container{T} = Union{AbstractVector{T}, NTuple{N,T}} where {N}

## ---- basic transformations needed for Mobius transform

# p -> p(-x)
function Base.reverse!(p::NTuple{N,T}) where {N,T}
    n = N ÷ 2
    for i ∈ 1:n
        j = N + 1 - i
        pᵢ, pⱼ = p[i], p[j]
        swap!(pᵢ, pⱼ)
    end
    p
end

# p -> p(λx)
function scale!(p::Container{T}, λ::T) where {T <: Union{BigFloat, 𝐈}}
    
    L = max(precision(first(p)), precision(λ))
    w, d, u = BigFloat(;precision=L),BigFloat(;precision=L),BigFloat(;precision=L)
    a = one(T)
    scale!(p, λ, a, w, d, u)
    
end

# non-allocating
function scale!(p::Container{T}, λ::T,
                a::T, w, d, u) where {T <: Union{BigFloat, 𝐈}}
    for i ∈ 2:length(p)
        pᵢ = p[i]
        rmul!(a, λ, w, d, u)
        rmul!(pᵢ, a, w, d, u)
    end
    
    return nothing
    
end

# used by lowerbound... which may see non-mutable Int or Float64 values
function scale!(p::Container{T}, λ::T)  where {T}
    a = one(T)
    n = length(p)
    for i ∈ 2:n
        a *= λ
        p[i] *= a

    end
    return nothing
end

# p -> p(x + λ)  This is 𝑶(n²) Horner scheme and could be faster
function taylor_shift!(ps::Container{T}, λ=one(T)) where {T<:𝐈}
    N = length(ps)
    L = maximum(precision, (ps[1].lo, ps[1].hi, λ.lo, λ.hi))
    # if we could reset precision, would be nice to not construct
    # this each time...
    w,d,u = BigFloat(;precision=L),BigFloat(;precision=L),BigFloat(;precision=L)
    taylor_shift!(ps, λ, ball(0, L), L, w, d, u)
end

# in place: pⱼ <- pⱼ + λ⋅pⱼ₊₁ sweeps; `tmp` is an interval used as storage
function taylor_shift!(ps::Container{T}, λ, tmp, L, w, d, u) where {T<:𝐈}
    N = length(ps)
    @inbounds for i in 1:N-1
        for j in N-1:-1:i
            copy!(tmp, ps[j+1])
            rmul!(tmp, λ, w, d, u)
            radd!(ps[j], tmp)
        end
    end
    return nothing
end

# p -> p(x + 1); only additions are needed
function taylor_shift1!(ps::Container{T}) where {T<:𝐈}
    N = length(ps)
    @inbounds for i in 1:N-1
        for j in N-1:-1:i
            radd!(ps[j], ps[j+1])
        end
    end
    return nothing
end


# p -> p((ax + b)/(x+b))
function mobius_transform!(p::Container{T}, a, b) where {T}
    
    L = max(precision(first(p)), precision(a), precision(b))
    w,d,u = BigFloat(;precision=L), BigFloat(;precision=L), BigFloat(;precision=L)
    tmp = ball(0, L)
    λ = ball(1, L)

    taylor_shift!(p, a, tmp, L, w, d, u)
    scale!(p, (b-a), λ, w, d, u)
    reverse!(p)
    taylor_shift1!(p)
    nothing
end

# basic structure of above without the storage setup to lessen allocations
function _mobius_transform!(p::Container{T}, a, b) where {T}
    taylor_shift!(p, a)
    scale!(p, (b-a))
    reverse!(p)
    taylor_shift1!(p)
    nothing
end      

## -----

# bound with interval arithmetic
function descartesbound(p, a::𝐈, b::𝐈, L::Int)
    ϵ = epsᴸ(L)
    q = ball.(deepcopy.(p), ϵ)
    descartesbound!(q, a, b, L)
end

function descartesbound!(q, a::𝐈, b::𝐈, L::Int)
    v, u, lead = descartesvar!(q, a, b, L)
    v + lead
end

# Sign variations of the Möbius transform of `q` for [a,b] at precision L.
# Returns (v, u, lead): `v` counts the variations among the coefficients
# whose sign is certain, `u` the number of coefficients whose sign is
# uncertain, and `lead` how many of those come before the first certain
# one. The true number of sign variations lies in [v, v + 2u].
function descartesvar!(q, a::𝐈, b::𝐈, L::Int)
    setprecision(L) do

        mobius_transform!(q, a, b)

        tol = epsᴸ(L) # 1/2^L
        tol_ = -tol
        i₀ = findfirst(q₀ -> q₀.hi <= tol_ || q₀.lo >= tol, q)
        i₀ == nothing && return (0, length(q), length(q))
        posflag = q[i₀].lo >= tol

        v = u = 0
        for i ∈ i₀:length(q)
            qᵢ = q[i]
            flag = qᵢ.lo >= tol ? true :
                qᵢ.hi <= tol_ ? false : nothing
            if flag == nothing # A "0"
                u += 1
                continue
            end
            if posflag != flag
                v += 1
                posflag = flag
            end
        end
        (v, u + i₀ - 1, i₀ - 1)
    end
end

# Descartes bounds for `p` on [a,m] and [m,b]. The working precision
# starts at L₀ and doubles until the signs of all the transformed
# coefficients are certain, or `done(va, vb)` reports that the
# variations already found settle the test. At the worst-case precision
# Lmax, coefficients of uncertain sign are treated as 0 (leading ones
# count as a variation). Also returns the last precision used.
function descartespair(p, a, m, b, L₀, Lmax, done)
    L = min(L₀, Lmax)
    while true
        ϵ = epsᴸ(L)
        𝐚, 𝐦, 𝐛 = ball.((a, m, b), ϵ)
        𝐪 = [ball(pᵢ, ϵ) for pᵢ ∈ p]

        va, ua, la = descartesvar!(𝐪, 𝐚, 𝐦, L)
        if L < Lmax && done(va, 0)
            return (va, 0, L)
        end
        ball!.(𝐪, p, ϵ)
        vb, ub, lb = descartesvar!(𝐪, 𝐦, 𝐛, L)

        L >= Lmax && return (va + la, vb + lb, L)
        ua == ub == 0 && return (va, vb, L)
        done(va, vb) && return (va, vb, L)
        L = min(2L, Lmax)
    end
end

