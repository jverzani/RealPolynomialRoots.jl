## lightweight IntervalArithmetic over BigFloat following paper
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

# 𝐱 <- [x↓ - ϵ, x↑ + ϵ]; the directed roundings make this an enclosure
# of x even when x is not representable at the precision of 𝐱
function ball!(𝐱::𝐈, x, ϵ::BigFloat)
    _setdown!(𝐱.lo, x)
    rsub!(𝐱.lo, ϵ, RoundDown)
    _setup!(𝐱.hi, x)
    radd!(𝐱.hi, ϵ, RoundUp)
    𝐱
end

_setdown!(y::BigFloat, x::BigFloat) = set!(y, x, RoundDown)
_setup!(y::BigFloat, x::BigFloat) = set!(y, x, RoundUp)
# exact conversions to BigFloat for these
_exactbig(x::Union{Float16, Float32, Float64}) = BigFloat(x, precision=53)
_exactbig(x::Integer) = BigFloat(x, precision=max(64, ndigits(x, base=2)))
for (f, dir) in ((:_setdown!, :prevfloat), (:_setup!, :nextfloat))
    @eval begin
        $f(y::BigFloat, x::Union{Integer, Float16, Float32, Float64}) = $f(y, _exactbig(x))
        # otherwise round to nearest, then step one ulp outward
        $f(y::BigFloat, x) = $f(y, $dir(BigFloat(x, precision=precision(y))))
    end
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


function swap!(a::T, b::T) where {T <: 𝐈}
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




    # case analysis on signs; only when both intervals contain 0
    # are more than 2 products needed
    if d₁ >= 0
        if d₂ >= 0
            mul!(d, d₁, d₂, RoundDown); mul!(u, u₁, u₂, RoundUp)
        elseif u₂ <= 0
            mul!(d, u₁, d₂, RoundDown); mul!(u, d₁, u₂, RoundUp)
        else
            mul!(d, u₁, d₂, RoundDown); mul!(u, u₁, u₂, RoundUp)
        end
    elseif u₁ <= 0
        if d₂ >= 0
            mul!(d, d₁, u₂, RoundDown); mul!(u, u₁, d₂, RoundUp)
        elseif u₂ <= 0
            mul!(d, u₁, u₂, RoundDown); mul!(u, d₁, d₂, RoundUp)
        else
            mul!(d, d₁, u₂, RoundDown); mul!(u, d₁, d₂, RoundUp)
        end
    else
        if d₂ >= 0
            mul!(d, d₁, u₂, RoundDown); mul!(u, u₁, u₂, RoundUp)
        elseif u₂ <= 0
            mul!(d, u₁, d₂, RoundDown); mul!(u, d₁, d₂, RoundUp)
        else
            mul!(d, d₁, u₂, RoundDown)
            mul!(working, u₁, d₂, RoundDown)
            min!(d, working)
            mul!(u, d₁, d₂, RoundUp)
            mul!(working, u₁, u₂, RoundUp)
            max!(u, working)
        end
    end

    set!(d₁, d, RoundDown)
    set!(u₁, u, RoundUp)

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
        return !(x.hi < -ϵ || x.lo > ϵ)
        x.lo < -ϵ && x.hi > ϵ
    end
end
