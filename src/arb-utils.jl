## our Intervals are of type ArbFloat{L} for possibly
## different Ls. This ensures promotion to largest L
_promote(x,y) = promote(x,y)
_promote(x::Int, y::Int) = (ArbFloat{53}(x),ArbFloat{53}(y))
_promote(x::Int, y::ArbFloat{L}) where {L} = (ArbFloat{L}(x),ArbFloat{L}(y))
_promote(x::ArbFloat{L}, y::Int) where {L} = (ArbFloat{L}(x),ArbFloat{L}(y))

function _promote(x::Float64, y::Float64)
    (ArbFloat{53}(x),ArbFloat{53}(y))
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

function _promote(x::T, y::BigFloat) where {L, T <: ArbFloat{L}}
    L′ = max(L, arbL(y))
    T′ = ArbFloat{L′}
    (T′(x), T′(y))
end

"same as ArbNumerics.workingprecision"
arbL(x::ArbFloat{L}) where {L} = L
arbL(x::BigFloat) = Base.MPFR.precision(x) #106
