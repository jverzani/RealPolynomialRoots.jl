# -----
## refinement

"""
    refine_roots(st::State)

Return estimates for the isolated roots using a bisection method due to Alefeld, Potra, and Shi.
"""
function refine_roots(st::State)
    λ = x -> evalpoly(x, st.p)
    [a42(λ, Tuple(I)) for I in st]
end

## From simple.jl of Roots
function a42(f, ab; atol=nothing, rtol=nothing, λ=0.7, μ=0.5)
    a, b = adjust_bracket(ab)
    δ₀ = b - a
    fa, fb = f(a), f(b)
    assert_bracket(fa, fb)

    tols = (
        λ   = λ,
        atol = isnothing(atol) ? zero(one(a)) : atol,
        rtol = isnothing(rtol) ? eps(one(a)) : rtol,
    )
    c = a - fa * (b - a) / (fb - fa)
    c = avoid_boundaries(a, c, b, fa, fb, tols)

    fc = f(c)
    iszero(fc) && return c
    e, fee = c, fc
    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    n = 2
    while true
        δ = tolₑ(a, b, fa, fb, tols.atol, tols.rtol)
        (b - a) ≤ δ && return (abs(fa) < abs(fb) ? a : b)

        ee, fee = d, fd
        for k in 1:2
            if n == 2 || iszero(_pairwise_prod(fa, fb, fd, fee))
                c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)
            else
                c = ipzero(a, b, d, ee, fa, fb, fd, fee)
                if (c <= a || b <= c)
                    c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)
                end
            end
            n += 1
            c = avoid_boundaries(a, c, b, fa, fb, tols)
            fc = f(c)
            iszero(fc) && return c

            ee, fee = d, fd
            a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

            δ = tolₑ(a, b, fa, fb, tols.atol, tols.rtol)
            (b - a) ≤ 2δ && return (abs(fa) < abs(fb) ? a : b)
        end

        n += 1

        u, fu = abs(fa) < abs(fb) ? (a, fa) : (b, fb)
        c = u - 2 * fu * (b - a) / (fb - fa)

        if 2abs(c - u) > (b - a)
            c = a / 2 + b / 2
        end

        c = avoid_boundaries(a, c, b, fa, fb, tols)
        fc = f(c)
        iszero(fc) && return c

        ee, fee = d, fd
        a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
        δ = tolₑ(a, b, fa, fb, tols.atol, tols.rtol)
        (b - a) ≤ 2δ && return (abs(fa) < abs(fb) ? a : b)

        if (b - a) ≥ μ * δ₀
            c = a / 2 + b / 2
            fc = f(c)
            iszero(fc) && return c
            ee, fee = d, fd
            a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
        end
        n += 1
    end
end


## helper function: floating point, sorted, finite
function adjust_bracket(x0)
    u, v = map(float, extrema(x0))
    if u > v
        u, v = v, u
    end
    isinf(u) && (u = nextfloat(u))
    isinf(v) && (v = prevfloat(v))
    u, v
end
@inline isbracket(fa, fb) = sign(fa) * sign(fb) < 0
assert_bracket(fx0, fx1) = isbracket(fx0, fx1) || throw(ArgumentError(bracketing_error))

## adjustment before calling bracket
function avoid_boundaries(a, c, b, fa, fb, tols)
    δ = tols.λ * tolₑ(a, b, fa, fb, tols.atol, tols.rtol)

    if (b - a) ≤ 4δ
        c = a / 2 + b / 2
    elseif c ≤ a + 2δ
        c = a + 2δ
    elseif c ≥ b - 2δ
        c = b - 2δ
    end
    c
end
## Brent-style tole from paper
function tolₑ(a, b, fa, fb, atol, rtol)
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    return 2 * u * rtol + atol
end

# assume fc != 0
## return a1,b1,d with a < a1 <  < b1 < b, d not there
@inline function bracket(a, b, c, fa, fb, fc)
    if isbracket(fa, fc)
        # switch b,c
        return (a, c, b, fa, fc, fb)
    else
        # switch a,c
        return (c, b, a, fc, fb, fa)
    end
end

# iterative quadratic solution to P(x) = 0 where P=f(a) + f[a,b]*(x-a) + f[a,b,d]*(x-a)*(x-b)
function newton_quadratic(a, b, d, fa, fb, fd, k::Int)
    A = _fabd(a, b, d, fa, fb, fd)
    B = _fab(a, b, fa, fb)

    (iszero(A) || !isfinite(A)) && return a - fa / B

    r = sign(A) * sign(fa) > 0 ? a : b

    for i in 1:k
        P = fa + B * (r - a) + A * (r - a) * (r - b)
        P′ = (B + A * (2r - a - b))
        r -= P / P′
    end

    return r
end

# f[a, b] divided differences
@inline _fab(a, b, fa, fb) = (fb - fa) / (b - a)

# f[a,b,d]
@inline function _fabd(a, b, d, fa, fb, fd)
    fab, fbd = _fab(a, b, fa, fb), _fab(b, d, fb, fd)
    (fbd - fab) / (d - a)
end

# check if fa,fb,fc,fd are distinct
function _pairwise_prod(as...)
    t = one(first(as))
    n = length(as)
    for i in 1:(n - 1)
        for j in (i + 1):n
            t *= (as[i] - as[j])
        end
    end
    t
end

# zero of inverse interpolation polynomial through (a,fa), (b,fb), (c,fc), (d, fd)
# may not lie in [a,b], though asymptotically will under smoothness assumptions
function ipzero(a, b, c, d, fa, fb, fc, fd)
    Q11 = (c - d) * fc / (fd - fc)
    Q21 = (b - c) * fb / (fc - fb)
    Q31 = (a - b) * fa / (fb - fa)
    D21 = (b - c) * fc / (fc - fb)
    D31 = (a - b) * fb / (fb - fa)
    Q22 = (D21 - Q11) * fb / (fd - fb)
    Q32 = (D31 - Q21) * fa / (fc - fa)
    D32 = (D31 - Q21) * fc / (fc - fa)
    Q33 = (D32 - Q22) * fa / (fd - fa)
    a + (Q31 + Q32 + Q33)
end
