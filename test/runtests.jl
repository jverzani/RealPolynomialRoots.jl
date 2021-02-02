using RealPolynomialRoots, Polynomials
using Test


@testset "Real roots" begin
    
    # produces a very close pair of roots
    function mignotte(n,τ)
        x = variable(SparsePolynomial)
        p = x^n - ((2^(τ ÷ 2)-1)*x - 1)^2
        convert(Polynomial, p)
    end

    st = real_roots(mignotte(7, 14))
    @test length(st.Isol) == 3
    @test length(st.Unresolved) == 0

    st = real_roots(mignotte(15, 14))
    @test length(st.Isol) == 3
    @test length(st.Unresolved) == 0


#    st = real_roots(mignotte(127, 14))# too slow
#    @test length(st.Isol) == 3
#    @test length(st.Unresolved) == 0

    p = mignotte(7, 14)
    x = variable(Polynomial{BigInt})
    st = real_roots(prod(p(x+i) for i in -3:3))
    @test length(st.Isol) == 7*3
    @test length(st.Unresolved) == 0


    q = prod(x^2 + i for i in 1:10)
    st = real_roots(p*q)
    @test length(st.Isol) == 3
    @test length(st.Unresolved) == 0

    ps = setprecision(1024) do
        delta = big(1e-100); 
        fromroots(Polynomial,[1.0, 1+delta, 1 + 2delta])
    end
    st = real_roots(ps)
    @test length(st.Isol) == 3
    @test length(st.Unresolved) == 0

    N = 50
    ps = [big(i)/N for i in 1:N]
    p = fromroots(Polynomial, ps);
    st = real_roots(p)
    @test length(st.Isol) == 50
    @test length(st.Unresolved) == 0


end
