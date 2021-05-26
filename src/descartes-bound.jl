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
    dps = [ball(0, L) for _ ∈ 1:N]
    taylor_shift!(ps, λ, dps, L, w, d, u)
end

function taylor_shift!(ps::Container{T}, λ, dps, L, w, d, u) where {T<:𝐈}

    N = length(ps)
    @inbounds for i in N:-1:1
        for j in N:-1:2
            j + i > N+1 && continue
            muladd!(dps[j], λ, dps[j-1], w, d, u)
        end
        muladd!(dps[1], λ, ps[i], w, d, u)
    end
    copy!.(ps, dps)
    return nothing
end


# p -> p((ax + b)/(x+b))
function mobius_transform!(p::Container{T}, a, b) where {T}
    
    L = max(precision(first(p)), precision(a), precision(b))
    w,d,u = BigFloat(;precision=L), BigFloat(;precision=L), BigFloat(;precision=L)
    dps = [ball(0, L) for _ ∈ 1:length(p)]
    λ = ball(1, L)
    
    taylor_shift!(p, a, dps, L, w, d, u);  zero!.(dps)
    scale!(p, (b-a), λ, w, d, u)
    reverse!(p)
    taylor_shift!(p, one(T), dps, L, w, d, u)
    nothing
end      

function Xmobius_transform!(p::Container{T}, a, b) where {T}
    taylor_shift!(p, a)
    scale!(p, (b-a))
    reverse!(p)
    taylor_shift!(p, one(T))
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
    setprecision(L) do

        mobius_transform!(q, a, b)

        tol = epsᴸ(L) # 1/2^L
        tol_ = -tol
        i₀ = findfirst(q₀ -> q₀.hi <= tol_ || q₀.lo >= tol, q)
        cnt = i₀-1 # count 0s
        q₀ = q[i₀]
        posflag = q₀.lo >= tol

        for i ∈ i₀:length(q)
            copy!(q₀, q[i])
            flag = q₀.lo >= tol ? true :
                q₀.hi <= tol_ ? false : nothing
            flag == nothing && continue # A "0"
            if posflag != flag
                cnt += 1
                posflag = flag
            end
        end
        cnt
    end
end

