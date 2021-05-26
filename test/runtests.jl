using RealPolynomialRoots
using Polynomials
using Test


@testset "Real roots" begin
    
    # produces a very close pair of roots
    # like 2^(-n*τ/2)
    function mignotte(n,τ)
        x = variable(SparsePolynomial)
        p = x^n - ((2^(τ ÷ 2)-1)*x - 1)^2
        coeffs(convert(Polynomial, p))
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
    
    𝐩 = Polynomial(mignotte(7, 14))
    𝐱 = variable(Polynomial{BigInt})

    
    p = coeffs(prod(𝐩(𝐱+i) for i in -3:3))
    st = ANewDsc(p)
    @test length(st) == 7*3
    @test length(st.Unresolved) == 0

    𝐪 = 𝐩 * prod(𝐱^2 + i for i in 1:10)
    p = coeffs(𝐪)
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0

    N = 50
    ps = [big(i)/N for i in 1:N]
    𝐩 = fromroots(Polynomial, ps);
    p = coeffs(𝐩)
    st = ANewDsc(p)
    @test length(st) == 50
    @test length(st.Unresolved) == 0

    # broken test
    # This passes with ^30
    𝐩 = setprecision(256) do
        delta = 1/big(10)^30
        fromroots(Polynomial,[1.0, 1+delta, 1 + 2delta])
    end
    p = coeffs(𝐩)
    st = ANewDsc(p)
    @test length(st) == 3
    @test length(st.Unresolved) == 0
    
    # This fails to disambiguate the roots, though
    # with ^50
    𝐩 = setprecision(256) do
        delta = 1/big(10)^50
        fromroots(Polynomial,[1.0, 1+delta, 1 + 2delta])
    end
    p = coeffs(𝐩)
    st = ANewDsc(p)
    @test_broken (length(st) == 3) && (length(st.Unresolved) == 0)

    # test real_roots
    rroots = (refine_roots ∘ ANewDsc ∘ coeffs)
    
    𝐱 = variable(Polynomial)

    # d = 0
    𝐩 = zero(𝐱)
    rts = rroots(𝐩)
    @test isempty(rts)

    𝐩 = one(𝐱)
    rts = rroots(𝐩)
    @test isempty(rts)


    # d = 1

    𝐩 = 𝐱
    rts = rroots(𝐩)
    @test length(rts) == 1 && abs(rts[1]) <= eps()

    𝐩 = 𝐱 + 1
    rts = rroots(𝐩)
    @test Float64.(rts) ≈ [-1]

    for i ∈ (1e2, 1e5, 1e8) # bigger values are just close
        𝐩 = 𝐱 - i
        rts = rroots(p)
        @test rts ≈ [i]
    end
    
    # d > 1, nreal = 1
    for i ∈ -2000:117:2000
        𝐩 = (𝐱^2+1)*(𝐱-i)
        rts = rroots(𝐩)
        @test rts ≈ [i]
    end

    # d = 2
    𝐩 = 𝐱*(𝐱-1)
    rts = rroots(𝐩)
    @test sort(Float64.(rts)) ≈ [0,1]


    # larger d
    xs = -10:2:10
    𝐩 = fromroots(Polynomial, xs)
    𝐱 = variable(𝐩)
    𝐪 = 𝐩 * (𝐱^2 + 1)*(𝐱^4 + 1)
    rts = rroots(𝐪)
    @test length(rts) == length(xs)
    @test maximum(abs, sort(rts) - xs) <= eps()

    


end
