using RealPolynomialRoots
using Test

# Minimal helpers for coefficient vectors [a₀, a₁, …, aₙ] (ascending order)
function polymul(p, q)
    T = promote_type(eltype(p), eltype(q))
    r = zeros(T, length(p) + length(q) - 1)
    for (i, a) in enumerate(p), (j, b) in enumerate(q)
        r[i+j-1] += a * b
    end
    r
end

function polyadd(p, q)
    T = promote_type(eltype(p), eltype(q))
    r = zeros(T, max(length(p), length(q)))
    r[1:length(p)] .+= p
    r[1:length(q)] .+= q
    r
end

# coefficients of (x - r₁)(x - r₂)⋯(x - rₙ)
fromroots(rs) = foldl((p, r) -> polymul(p, [-r, one(r)]), rs; init=[one(eltype(rs))])

# coefficients of p(x + c), by Horner's method
function pshift(p, c)
    r = [p[end]]
    for a in reverse(p[1:end-1])
        r = polyadd(polymul(r, [c, one(c)]), [a])
    end
    r
end

@testset "Real roots" begin

    # produces a very close pair of roots
    # like 2^(-n*τ/2)
    # x^n - ((2^(τ÷2) - 1)x - 1)^2
    function mignotte(n,τ)
        a = 2^(τ ÷ 2) - 1
        p = zeros(Int, n + 1)
        p[end] = 1
        p[1:3] .-= [1, -2a, a^2]
        p
    end

    p = mignotte(7, 14)
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0

    p = mignotte(15, 14)
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0

    p = mignotte(63, 14)    # gets slow ~ 2secs
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0

    𝐩 = mignotte(7, 14)

    p = reduce(polymul, (pshift(big.(𝐩), big(i)) for i in -3:3))
    st = ANewDsc(p)
    @test length(st) == 7*3
    @test length(st.Unresolved) == 0

    p = reduce(polymul, ([i, 0, 1] for i in 1:10); init=𝐩)
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0

    N = 50
    ps = [big(i)/N for i in 1:N]
    p = fromroots(ps)
    st = ANewDsc(p)
    @test length(st) == 50
    @test length(st.Unresolved) == 0

    # These two tests show that polynomials may not have the roots
    # that they theoretically might expect to have due to round
    # off.
    # This passes with ^30
    p = setprecision(256) do
        delta = 1/big(10)^30
        fromroots([1.0, 1+delta, 1 + 2delta])
    end
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0

    # This shows mathematical expectations may fail due to precision issues
    delta = 1e-5
    p = fromroots([1.0, 1+delta, 1 + 2delta])
    st = ANewDsc(p)
    @test length(st) == 3 # works

    delta = 1e-8
    p = fromroots([1.0, 1+delta, 1 + 2delta])
    st = ANewDsc(p)
    @test_broken length(st) == 3 # doesn't to find 3

    delta = 1e-8
    p = fromroots(BigFloat[1.0, 1+delta, 1 + 2delta])
    st = ANewDsc(p)
    @test length(st) == 3 # works

    ## test real_roots
    rroots = (refine_roots ∘ ANewDsc)

    # d = 0
    p = [0.0]
    rts = rroots(p)
    @test isempty(rts)

    p = [1.0]
    rts = rroots(p)
    @test isempty(rts)


    # d = 1

    p = [0.0, 1.0]
    rts = rroots(p)
    @test length(rts) == 1 && abs(rts[1]) <= eps()

    p = [1.0, 1.0]
    rts = rroots(p)
    @test Float64.(rts) ≈ [-1]

    for i ∈ (1e2, 1e5, 1e8) # bigger values are just close
        p = [-i, 1.0]
        rts = rroots(p)
        @test Float64.(rts) ≈ [i] atol=1
    end

    # d > 1, nreal = 1
    for i ∈ -2000:117:2000
        p = polymul([1.0, 0.0, 1.0], [-i, 1.0])
        rts = rroots(p)
        @test Float64.(rts) ≈ [i] atol=abs(i)
    end

    # d = 2
    p = [0.0, -1.0, 1.0]
    rts = rroots(p)
    @test sort(Float64.(rts)) ≈ [0,1]


    # larger d
    xs = -10:2:10
    p = fromroots(collect(xs))
    q = polymul(polymul(p, [1, 0, 1]), [1, 0, 0, 0, 1])
    rts = rroots(q)
    @test length(rts) == length(xs)
    @test maximum(abs, sort(rts) - xs) <= eps()

end
