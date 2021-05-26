
## ----
##
## lightweight IntervalArithemetic over BigFloat following paper
## (https://arxiv.org/pdf/1104.1362.pdf) on root refinement
## https://github.com/JuliaIntervals/IntervalArithmetic.jl is the
## proper source for such calculations
struct 𝐈{T}
    lo::T
    hi::T
end

## An interval can be iterated over to return the end points
Base.length(I::𝐈) = 2
function Base.iterate(I::𝐈, state=nothing)
    isnothing(state) && return (I.lo, 1)
    isone(state) && return (I.hi, 2)
    return nothing
end        

function Base.show(io::IO, I::𝐈)
    a, b = I
    L = maximum(precision, (a,b))
    w = b-a

    w <= 0 && (print(io, "{$a}"); return nothing)
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
    nothing
end


𝑰(a,b,args...) = error("");#𝐈(a,b)

Base.diff(a::𝐈) = a.hi-a.lo
midpoint(a::𝐈) = a.lo + (a.hi - a.lo)/2
Base.sign(a::𝐈) = a.hi < 0 ? -1 : a.lo > 0 ? 1 : 0
minabs(a::𝐈) = a.hi < 0 ? -a.hi : (a.lo > 0) ? a.lo : zero(a.lo)
maxabs(a::𝐈) = max(abs(a.lo), abs(a.hi))
Base.isinf(a::𝐈) = isinf(a.lo) || isinf(a.hi)
Base.isnan(a::𝐈) = isnan(a.lo) || isnan(a.hi)
function Base.zero(a::𝐈)
    L = maximum(precision, (a.lo, a.hi))
    ball(0,L)
end
zero!(a::𝐈) = (zero!(a.lo); zero!(a.hi); nothing)
    

precision(I::𝐈) = maximum(precision, I)


function ball!(𝐱::𝐈, x, ϵ::BigFloat)
    x′ = BigFloat(x, precision=precision(ϵ)) #deepcopy(x)
    set!(𝐱.lo, x′, RoundDown)
    rsub!(𝐱.lo, ϵ, RoundDown)
    set!(𝐱.hi, x′, RoundUp)
    radd!(𝐱.hi, ϵ, RoundUp)
    𝐱
end

function ball!(𝐱, x, L::Int)
    ϵ = epsᴸ(L)
    ball!(𝐱, x, ϵ)
end

function ball(x, ϵ::BigFloat)
    𝐱 = ball(precision(ϵ))
    ball!(𝐱, x, ϵ)
end

function ball(x,L::Int)
    𝐱 = ball(L)
    ball!(𝐱, x, L)
end

function ball(x::Int,L::Int)
    T = BigFloat
    return 𝐈(T(x, precision=L), T(x, precision=L))
    𝐱 = ball(L)
    x′ = T(x, precision=L) #deepcopy(x)
    ball!(𝐱, x′, L)
    𝐱
end

ball(L::Int) = 𝐈(BigFloat(;precision=L), BigFloat(;precision=L))
ball(ϵ::BigFloat) = I(BigFloat(;precision=precision(ϵ)), BigFloat(;precision=precision(ϵ)))

Base.one(::Type{T}) where {T <: 𝐈} = 𝐈(BigFloat("1"), BigFloat("1"))

function Base.copy!(a::𝐈, b::𝐈)
    set!(a.lo, b.lo)
    set!(a.hi, b.hi)
end





function swap!(a::𝐈, b::𝐈)
    swap!(a.lo, b.lo)
    swap!(a.hi, b.hi)
end

function Base.:+(a::𝐈, b::𝐈)
    c = deepcopy(a)
    radd!(c, b)
    c
end

function radd!(a::𝐈, b::𝐈)
    d₁,u₁ = a.lo, a.hi
    d₂,u₂ = b.lo, b.hi
    T = BigFloat
    radd!(d₁, d₂, RoundDown)
    radd!(u₁, u₂, RoundUp)    
end

function radd!(a::𝐈, b::Number)
    T = BigFloat
    radd!(a.lo, b, RoundDown)
    radd!(a.hi, b, RoundUp)
end

function Base.:-(a::𝐈)
    d₁,u₁ = a.lo, a.hi
    𝐈(-a.hi, -a.lo)
end

Base.:-(a::𝐈, b::𝐈) = a + (-b)

function p(a::𝐈, b::𝐈)
    down = d1d2, u1u2, u1d2, d1u2
    up   = d1d2, u1u2, u1d2, d1u2

    d2d1
end


# right multiplication a <- a*b; 
# needs three BigFloat numbers to be non-allocating
# working, d, u make non-allocating
function LinearAlgebra.rmul!(a::𝐈, b::𝐈, working′=nothing, d′=nothing, u′=nothing) # working a bigplot
    d₁,u₁ = a.lo, a.hi
    d₂,u₂ = b.lo, b.hi

    L = maximum(precision, (d₁,u₁,d₂,u₂))
    working = working′ != nothing ? working′ : BigFloat(;precision=L)
    d = d′ == nothing ? BigFloat(;precision=L) : d′
    u = u′ == nothing ? BigFloat(;precision=L) : u′



    
    set!(d, d₁); set!(u, d₁)
    rmul!(d, d₂, RoundDown); rmul!(u, d₂, RoundUp)
    
    
    set!(working, u₁) 
    rmul!(working, u₂, RoundDown)
    min!(d, working)

    set!(working, u₁)
    rmul!(working, d₂, RoundDown)
    min!(d, working)

    set!(working, d₁)
    rmul!(working, u₂, RoundDown)
    min!(d, working)

    set!(working, d₁)
    rmul!(working, u₂, RoundUp)
    max!(u, working)

    set!(working, u₁)
    rmul!(working, d₂, RoundUp)
    max!(u, working)

    set!(working, u₁)
    rmul!(working, u₂, RoundUp)
    max!(u, working)

    set!(d₁,d)
    set!(u₁,u)

    a
end

function LinearAlgebra.rmul!(a::𝐈, b::Number)
    T = BigFloat
    rmul!(a.lo, b, RoundDown)
    rmul!(a.hi, b, RoundUp)
    
    if b < 0
        tmp = a.lo
        set!(a.lo, a.hi)
        set!(a.hi, tmp)
    end
    
    a
end

function inv!(b::𝐈)
    T = BigFloat
    c = deepcopy(b.lo)
    L = precision(c)
    set!(b.lo, rdiv!(T(1,precision=L), b.hi, RoundDown))
    set!(b.hi, rdiv!(T(1,precision=L), c, RoundUp))
    b
end
    
function LinearAlgebra.rdiv!(a::𝐈, b::𝐈)
    c = deepcopy(b)
    inv!(c)
    rmul!(a, c)
end

function swap!(c::T, a::T,b::T) where {T <: 𝐈}
    swap!(a.lo, b.lo)
    swap!(a.hi, b.hi)
end

function abs!(a::𝐈)
    if a.hi < 0
        c = deepcopy(a.lo)
        swap!(a.lo, -a.hi)
        swap!(a.hi, -c)
    end
end

# muladd with extra storage to pass to `rmul!`
function muladd!(a::𝐈, b::𝐈, c::𝐈, w=BigFloat(), d=BigFloat(), u=BigFloat())
    rmul!(a, b, w, d, u)
    radd!(a, c)
    nothing
end

# 0 is uncertain
maybe0(x::𝐈) = x.lo < zero(x.lo) < x.hi
function maybe0(x::𝐈, L::Int)
    setprecision(L) do
        ϵ = eps(BigFloat)
#        return x.hi > -ϵ && x.lo < ϵ
        return !(x.hi < -ϵ || x.lo > ϵ)
        x.lo < -ϵ && x.hi > ϵ
    end
end
