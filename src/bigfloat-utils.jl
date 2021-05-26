## -----
## BigFloat numbers are used; The MutableArithmetics.jl package was
## mimicked below, though not included directly, so rounding modes
## could be specified.  This implementation below still allocates much
## more than is desirable.


# for more generic uses
Big(a;precision=53) = BigFloat(a; precision=precision)

#import Base: precision
precision(x::Float64) = 53
precision(x::BigFloat) = Base.MPFR.precision(x)

## ---- bigfloat
## We add some methods from libmpfr to give rounding
## hints. (cf. MutableArithmetics.jl for inspiration)

# This caused errors which is too bad
# as it means many more allocations because we must create new
# BigFloat values every time the precision is changes
set_prec!(x::BigFloat,n::Int) = ccall(
    (:mpfr_set_prec, :libmpfr),
    Cvoid,
    (Ref{BigFloat}, Cint, Base.MPFR.MPFRRoundingMode),
    x,Cint(n), Base.MPFR.ROUNDING_MODE[]
)

# this is copy!(y, x)
# we can not use `set_prec` to set precision though
set!(y::BigFloat, x, rnd=Base.MPFR.ROUNDING_MODE[]) = ccall(
    (:mpfr_set, :libmpfr),
    Int32,
    (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
    y, x, rnd
)
zero!(a::BigFloat) = set!(a, 0)

swap!(y::BigFloat,x::BigFloat) = ccall(
    (:mpfr_swap, :libmpfr),
    Cvoid,
    (Ref{BigFloat}, Ref{BigFloat}),
    y, x
)

# in place +,-,*./
function radd!(a::BigFloat, b::BigFloat, rnd=Base.MPFR.ROUNDING_MODE[])
    ccall(
        (:mpfr_add, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        a,a,b, rnd
    )
    a
end

function rsub!(a::BigFloat, b::BigFloat, rnd=Base.MPFR.ROUNDING_MODE[])
    ccall(
        (:mpfr_sub, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        a,a,b, rnd
    )
    a
end

function LinearAlgebra.rmul!(a::BigFloat, b::Number, rnd=Base.MPFR.ROUNDING_MODE[])
    ccall(
        (:mpfr_mul, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        a,a,b, rnd
    )
    a
end

function LinearAlgebra.rdiv!(a::BigFloat, b::BigFloat, rnd=Base.MPFR.ROUNDING_MODE[])
    ccall(
        (:mpfr_div, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        a,a,b, rnd
    )
    a
end

# a <- a*b+c
function muladd!(a::T,b,c) where {T <: BigFloat}
    rmul!(a,b)
    radd!(a,c)
end

# nonallocating online minimum and maximum
function min!(x::T,a::T) where {T <: BigFloat}
    a < x && set!(x, a)
    nothing
end

function max!(x::T,a::T) where {T <: BigFloat}
    a > x && set!(x, a)
    nothing
end

epsᴸ(L::Int) = setprecision(() -> eps(BigFloat), L)

# evalpoly with interval arithmetic to precision L
function evalpolyᴸ(x, p, L)
    Σ = ball(0, L)    
    evalpolyᴸ!(Σ, x, p, L)
    return Σ
end

function evalpolyᴸ!(Σ, x, p, L)
    ϵ = epsᴸ(L)
    𝐱 = ball(x, ϵ)
    𝐩ᵢ = ball(L)

    w,d,u = BigFloat(;precision=L),BigFloat(;precision=L),BigFloat(;precision=L)

    evalpolyᴸ!(Σ, 𝐱, p, L, 𝐩ᵢ, w, d, u)
end

function evalpolyᴸ!(Σ, 𝐱,  p, L, # storage, fat x, p, L
                    𝐩ᵢ, w, d, u)  # storage

    ϵ = epsᴸ(L)
    ball!(𝐩ᵢ, p[end], ϵ) 
    setprecision(L) do
        copy!(Σ, 𝐩ᵢ)
        for i ∈ length(p)-1:-1:1
            # Σ = muladd(Σ, x, p[i])
            rmul!(Σ, 𝐱, w, d, u)
            ball!(𝐩ᵢ, p[i], ϵ)
            radd!(Σ, 𝐩ᵢ) 
        end
        Σ
    end
end
