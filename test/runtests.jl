using RealPolynomialRoots, Polynomials
using Test


@testset "Real roots" begin
    
    # produces a very close pair of roots
    function mignotte(n,τ)
        x = variable(SparsePolynomial)
        p = x^n - ((2^(τ ÷ 2)-1)*x - 1)^2
        convert(Polynomial, p)
    end

    p = mignotte(7, 14)    
    rts = real_roots(p)
    @test length(rts) == 3
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test length(st.Unresolved) == 0

    p = mignotte(15, 14)    
    rts = real_roots(p)
    @test length(rts) == 3
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test length(st.Unresolved) == 0

    p = mignotte(63, 14)    
    rts = real_roots(p)
    @test length(rts) == 3
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test length(st.Unresolved) == 0
    

#    st = real_roots(mignotte(127, 14))# too slow
#    @test length(st.Isol) == 3
#    @test length(st.Unresolved) == 0

    x = variable(Polynomial{BigInt})

    p = mignotte(7, 14)
    q = prod(p(x+i) for i in -3:3)
    rts = real_roots(q)
    @test length(rts) == 7*3
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test length(st.Unresolved) == 0


    q = prod(x^2 + i for i in 1:10)
    rts = real_roots(p*q)
    @test length(rts) == 3
    st = RealPolynomialRoots.ANewDsc(coeffs(p*q))
    @test length(st.Unresolved) == 0

    N = 50
    ps = [big(i)/N for i in 1:N]
    p = fromroots(Polynomial, ps);
    rts = real_roots(p)
    @test length(rts) == 50
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test length(st.Unresolved) == 0


end


@testset "Broken tests" begin

    N = 20
    p = setprecision(1024) do
        delta = (big(1)/10)^N
        fromroots(Polynomial,[1.0, 1+delta, 1 + 2delta])
    end
    rts = real_roots(p)
    @test length(rts) == 3
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test length(st.Unresolved) == 0

    ## Careful, N=30 gets into a linear mode trap
    
    ##XXX THIS IS FAILING XXX
    ##XXX Should be disambiguated, isn't XXX
    N = 100
    p = setprecision(1024) do
        delta = (big(1)/10)^100
        fromroots(Polynomial,[1.0, 1+delta, 1 + 2delta])
    end
    rts = real_roots(p)
    @test_broken length(rts) == 3
    st = RealPolynomialRoots.ANewDsc(coeffs(p))
    @test_broken length(st.Unresolved) == 0

end
