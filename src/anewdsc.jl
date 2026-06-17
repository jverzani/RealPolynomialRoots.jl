## ----
# Find isolating intervals for real roots of a polynomial
# following
# Computing real roots of real polynomials
# Michael Sagraloff, Kurt Mehlhorn
# https://doi.org/10.1016/j.jsc.2015.03.004




## ---- Admissible point

# pseudo admissible point uses randomization to cut down number of operations
#
# This simple function allocates too much:
# * the evalpoly call allocates
# * the allocations in evalpoly are due to allocations in the interval
#   arithmetic implementation of `rmul!`.
# * these allocations can be removed by passing in 3 temporary variables
# * these can't be recycled here, as we don't have a mechanism to change
#   the working precision of a BigFloat number. (`set_prec!` will cause
#   segmentation faults)
# This is an issue then for: evalpoly and the various tests
# zerotest, onetest, newtontest, boundarytest
function admissiblepoint(p, m, ϵ, n, L′ = 53)
    λ = 2 + ceil(Int, log2(n))
    L = L′
    L = max(L′, maximum(precision∘float, p), precision(m))
    while true
        n′ = λ*n
        i = rand(-n′:n′)
        δ = epsᴸ(L-2) # 4*2⁻ᴸ
        ξ = setprecision(() -> m + i * ϵ/n′, L)
        𝐯ᵢ = evalpolyᴸ(ξ, p, L)
        (𝐯ᵢ.lo > δ || 𝐯ᵢ.hi < -δ) && return ξ
        L = 2L
    end

end


## ---- Tests (onetest, zerotest, newtontest, boundarytest, signtest)

# log2(a)
function tₐ(pa,L=precision(pa))
    apa = abs(pa)
    ϵ = epsᴸ(L)
    apa > ϵ ? ceil(Int, log2(apa)) : 2L
end

τ(p::Container{T}) where {T} = max(one(T), log2(norm(p, Inf)))
Mlog(x::T) where {T} = log2(max(one(T), abs(x)))




## Test if the interval 𝑰=[a,b] has a zero in it
function zerotest(p, I)
    isempty(p) && return false
    a,b = I
    a >= b && return true

    L′ = max(maximum(precision∘float, p), maximum(precision, (a,b)))

    ta, tb = tₐ(a,L′), tₐ(b,L′)

    n = length(p) - 1
    n′ = sum(divrem(n,2))

    L′′ = max(24, L′, max(1, -min(ta-1, tb-1) + 2(n+1) + 1))
    L = ceil(Int, n + τ(p) + n * Mlog(a) + n * Mlog(b) + L′′)

    m = setprecision(() -> a + (b-a)/2, L)
    ϵ = epsᴸ(L)
    𝐚, 𝐦, 𝐛 = ball.((a,m,b), ϵ)

    𝐪 = ball.(deepcopy.(p), ϵ)
    bnda = descartesbound!(𝐪, 𝐚, 𝐦, L)
    !iszero(bnda) && return false


    ball!.(𝐪, p, ϵ)
    bndb = descartesbound!(𝐪, 𝐦, 𝐛, L)
    iszero(bndb) && return true

    return false

end

## test if the interval 𝑰 = [a,b] has exactly 1 root in it
function onetest(p, I)
    isempty(p) && return false
    a,b = I
    a >= b && return false

    L′ = max(maximum(precision∘float, p), maximum(precision, (a,b)))

    n = length(p) - 1
    n′ = sum(divrem(n,2))

    ϵ = (b-a) * 2.0^(-ceil(Int, log2(n)+2)) #   (b-a) / (4n)
    m = setprecision(() ->a + (b-a)/2, L′)

    mstar = admissiblepoint(p, m, ϵ, n′)

    ta, tb,t = tₐ(a,L′), tₐ(b,L′),tₐ(mstar,L′)
    L′′ = max(24, L′, max(1, -min(ta-1,tb-1, t-1)) + 4(n+2))
    L = ceil(Int, n + τ(p) + n * Mlog(a) + n * Mlog(b) + L′′)

    ϵ = epsᴸ(L)
    𝐚,𝐦,𝐛 = ball.((a ,mstar, b), ϵ)
    𝐪 = ball.(deepcopy.(p), ϵ)

    bnda = descartesbound!(𝐪, 𝐚, 𝐦, L)

    ball!.(𝐪, p, ϵ)
    bndb = descartesbound!(𝐪, 𝐦, 𝐛, L)

    (bnda == -1 || bndb == -1) && throw(ArgumentError("L not big enough for descartesbound"))
    iszero(bnda) && isone(bndb) && return true
    isone(bnda) && iszero(bndb) && return true

    return false

end


# The newton test over I=[a,b] succeeds (informally) if:
# * a cluster of k roots (k ≈ λ below) is in a smaller interval J=[c,d]
# * there are no other roots nearby J
# Success means a new interval with size ≈≤ w(I)/N is returned
# Lemma 20 has this condition on success
# w(J) ≤ 2⁻¹³ w(I)/N
# k roots are in the one-circe region of Δ(J)
# the disk centered at m(I) with radius 2^(log(n)+10) * N * w(I) contains no other roots
function newtontest(p, p′, I, nᴵ)

    a, b = I
    nᴵ = min(5, nᴵ)
    N = 2^(2^nᴵ)

    n = length(p) - 1
    n′ = sum(divrem(n, 2))
    ϵ = 2.0^(-ceil(Int, (5 + log(n))))
    L::Int = Lᵐ = maximum(precision, (a,b))


    for (i,j) ∈ ((1,2),(1,3),(2,3))
        discard_ij = false
        double_cnt = 0
        L = Lᵐ
        while !discard_ij

            # double L until either (25) or (26) is satisfied for the pair i,j
            m, w = setprecision(() -> (a + (b-a)/2, b-a), L)
            uᵢ,uⱼ = setprecision(() -> (a + i*w/4,a + j*w/4), L)
            ξᵢ = admissiblepoint(p, uᵢ, ϵ * w, n′, L)
            ξⱼ = admissiblepoint(p, uⱼ, ϵ * w, n′, L)

            𝐩ᵢ, 𝐩ⱼ  = evalpolyᴸ(ξᵢ, p, L),  evalpolyᴸ(ξⱼ, p, L)
            𝐩′ᵢ,𝐩′ⱼ = evalpolyᴸ(ξᵢ, p′, L), evalpolyᴸ(ξⱼ, p′, L)

            L′ = max(L, maximum(precision, (ξᵢ, ξⱼ)))

            rdiv!(𝐩ᵢ,𝐩′ᵢ); rdiv!(𝐩ⱼ, 𝐩′ⱼ)
            𝚫 = 𝐩ᵢ - 𝐩ⱼ
            𝐯ᵢ, 𝐯ⱼ =  𝐩ᵢ, 𝐩ⱼ

            pᵢ, pⱼ = maxabs(𝐩ᵢ), maxabs(𝐩ⱼ) # store, used to compute new L

            if (min(minabs(𝐯ᵢ), minabs(𝐯ⱼ)) > w) || (4n*maxabs(𝚫) < w)  # (25)
                # discard_too_big_difference_too_small
                discard_ij = true # stop checking this pair

            elseif (max(maxabs(𝐯ᵢ), maxabs(𝐯ⱼ)) < 2w) && (8n*minabs(𝚫) > w) # (26)

                discard_ij = true # stop checking this pair

                # now set L to be big enough and estimate λ
                # L = ....²¹ (leaves diff between L an dactual < 1/(32N)

                L = max(24, ceil(Int, log2(n) + log2(N) +
                          Mlog(1/pᵢ) + Mlog(1/pⱼ) + max(0, -log2(w))))

                𝐩ᵢ,  𝐩ⱼ  = evalpolyᴸ(ξᵢ, p, L),  evalpolyᴸ(ξⱼ, p, L)
                𝐩′ᵢ, 𝐩′ⱼ = evalpolyᴸ(ξᵢ, p′, L), evalpolyᴸ(ξⱼ, p′, L)


                rdiv!(𝐩ᵢ,𝐩′ᵢ); rdiv!(𝐩ⱼ,𝐩′ⱼ)
                𝐯ᵢ,𝐯ⱼ = 𝐩ᵢ,𝐩ⱼ

                λ = 𝐯ᵢ
                rdiv!(λ, 𝐯ᵢ - 𝐯ⱼ)
                rmul!(λ, ξⱼ - ξᵢ)
                radd!(λ, ξᵢ)


                if λ.hi < a || λ.lo > b
                    continue
                elseif isinf(λ) || isnan(λ)
                    continue
                else

                    lᵢⱼ = setprecision(() ->floor(Int, (midpoint(λ) - a)/w * 4N), L)
                    aᵢⱼ = setprecision(() -> a + (lᵢⱼ - 1)  * w/(4N), L)
                    bᵢⱼ = setprecision(() -> a + (lᵢⱼ + 2) * w/(4N), L)

                    aᵢⱼᵅ = admissiblepoint(p, aᵢⱼ, ϵ *w/N, n′, L)
                    bᵢⱼᵅ = admissiblepoint(p, bᵢⱼ, ϵ *w/N, n′, L)

                    lᵢⱼ - 1 <= 0  && (aᵢⱼᵅ = Big(a, precision=L))
                    lᵢⱼ + 2 >= 4N && (bᵢⱼᵅ = Big(b, precision=L))

                    aᵢⱼᵅ, bᵢⱼᵅ = aᵢⱼᵅ < bᵢⱼᵅ  ? (aᵢⱼᵅ, bᵢⱼᵅ) : (bᵢⱼᵅ, aᵢⱼᵅ)

                    b-a < bᵢⱼᵅ - aᵢⱼᵅ && continue

                    bnda = zerotest(p, (a,  aᵢⱼᵅ ))
                    bndb = zerotest(p, (bᵢⱼᵅ, b))

                    if  bnda && bndb
                        return true, 𝐈(aᵢⱼᵅ,  bᵢⱼᵅ)
                    end

                end
            elseif double_cnt > 4
                # abandon and use  bisection
                discard_ij = true
            end

            double_cnt += 1
            L = 2L
        end
    end
    return false, 𝐈(a,b)

end

# The boundary test tests for a cluster of roots near an endpoint of I
# that is not well separated from other roots in Δ(I)
#
# return (true, J) can reduce interval with J ≈≤ w(I)/N
# return (false, I) can't trim, use linear
# no root in I should throw an error
function boundarytest(p, I, nᴵ)
    a, b = I
    nᴵ = min(5, nᴵ)
    N = 2^(2^nᴵ)

    n = length(p) - 1
    n′ = sum(divrem(n,2))
    ϵ = 1/2.0^(2 + ceil(Int, log(n)))

    L = maximum(precision, (a,b))

    w, m, mₗ, mᵣ = setprecision(L) do
        w = b-a
        m = a + w/2
        mₗ, mᵣ = a + w/(2N), b - w/(2N)
        w, m, mₗ, m
    end

    ϵ′ = ϵ*w/N
    mₗ⁺, mᵣ⁺ = admissiblepoint(p, mₗ, ϵ′, n′), admissiblepoint(p, mᵣ, ϵ′, n′)

    zₗ, zᵣ = zerotest(p, (mₗ⁺, b)), zerotest(p, (a, mᵣ⁺))
    zₗ && zᵣ && return (false, 𝐈(a,b)) #error("No root in I $a, $mᵣ⁺, $mₗ⁺  $b")
    zₗ && return (true, 𝐈(a, mₗ⁺))
    zᵣ && return (true, 𝐈(mᵣ⁺, b))
    return (false, 𝐈(a,b))

end

# true if `p` crosses in interval `I`
# false if not
# nothing if ambiguous
function signtest(p, I)
    a,b = I
    L = maximum(precision, (a,b))

    Σ = ball(0, L)
    𝐚 = evalpolyᴸ!(Σ, a, p, L)
    sa = 𝐚.hi < 0 ? -1 : 𝐚.lo > 0 ? 1 : return nothing
    𝐛 = evalpolyᴸ!(Σ, b, p, L)
    sb = 𝐛.hi < 0 ? -1 : 𝐛.lo > 0 ? 1 : return nothing

    sa * sb < 0

end

# split interval at an admissible point
function  bisect_interval(p, I, n′, L)
    setprecision(L) do
        a,b = I
        c = setprecision(L) do
            w = b-a
            m = 4 + ceil(Int, log2(n′))
            admissiblepoint(p, a + w/2, w/2^m, n′,L)
        end

        𝐈(a,c), 𝐈(c,b)
    end
end


## ----

## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
## δ is buffer ensure interval endpoint do not include roots
function upperbound(p, δ = 1/2)

    T = BigFloat
    n = length(p) - 1
    n > 0 || return δ
    L = 2 + n + ceil(Int, τ(p)) #+ 53

    setprecision(L) do
        descartescount(p) == 0 && return zero(T) + δ
        q, d = p./p[end], n

        d == 0 && return δ
        d == 1 && return max(0, -q[1]) + δ

        a1 = abs(q[d])
        B = maximum([abs(q[i]) for i in 1:(d-1)])

        a,b,c = 1, -(1+a1), a1-B
        out = (-b + sqrt(b^2 - 4a*c))/2
        out + δ
    end

end

function lowerbound(p)
    q = [(isodd(i) ? 1 : -1)*p[i] for i ∈ eachindex(p)]
    return -upperbound(q)
end

root_bound(p) = (lowerbound(p), upperbound(p))


# Descartes bound on real roots in [0,oo) based on the sign changes
# treat 0 (|aᵢ| <= eps(T)) as no sign change
function descartescount(q::Container{T}, tol::T=zero(T)) where {T}
    n = length(q)

    i₀ = 1 # leading 0s
    while abs(q[i₀]) <= tol
        i₀ += 1
    end

    q₀ = q[i₀]
    flag = q₀ < 0 ? true : false

    cnt = 0

    for i ∈ (i₀+1):n
        qᵢ = q[i]
        abs(qᵢ) <= tol && continue
        flag′ = qᵢ < 0 ? true : false
        if flag′ != flag
            cnt += 1
            flag = flag′
        end
    end

    cnt

end

## ----

# hold the state
struct State{N,T,S}
    Isol::Vector{𝐈{S}}                        # DesBound == 1
    Unresolved::Vector{𝐈{S}}
    p::NTuple{N,T}
end

(st::State)(x) = evalpoly(x, st.p)

## iterate over Isol
Base.length(st::State) = length(st.Isol)
Base.iterate(st::State, state=nothing) = state==nothing ? iterate(st.Isol) : iterate(st.Isol, state)

function Base.show(io::IO, st::State)

    if !isempty(st.Unresolved)
        println(io, "There are unresolved intervals:")
        println.(Ref(io), st.Unresolved)
        println(io, "")
    end

    n = length(st.Isol)

    msg = n == 0 ? "No isolating intervals found." :
        n == 1 ?  "There was 1 isolating interval found:" :
        "There were $n isolating intervals found:"

    println(io, msg)
    n > 0 && println.(Ref(io), st.Isol)

    return nothing
end

## ----

"""
    ANewDsc(p; root_bound=(lowerbound(p), upperbound(p)), max_depth=96)
    refine_roots(st::State)

A method to find isolating intervals for the real roots of a
square-free polynomial specified by the coefficients stored in `p`.

* `p`: the polynomial coefficients, `[a₀, a₁, …, aₙ]`, of a **square-free** polynomial.
* `root_bound`: a lower bound and upper bound for the real roots.


Returns a `State` instance which has components:

* `Isol` holding the isolating intervals. Iteration over a `State` object will iterate over `Isol`.
* `Unresolved` holding any unresolved intervals. The show method alerts the presence of any such intervals.

The algorithm has a random step included, which leads to small variations in the output.

Examples:

```
julia> using RealPolynomialRoots

julia> ps = [-1, 254, -16129, 0, 0, 0, 0, 1] # mignotte polynomial with two nearby roots
8-element Vector{Int64}:
     -1
    254
 -16129
      0
      0
      0
      0
      1

julia> st = ANewDsc(ps)
There were 3 isolating intervals found:
[3.5…, 7.25…]₅₃
[0.00787401594698…, 0.00787401779189…]₆₃
[0.00787401419348…, 0.00787401594698…]₆₃

julia> ps = [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]; # πᵢ₌₁¹⁰ (x-i)

julia> st = ANewDsc(ps)
There were 10 isolating intervals found:
[9.25…, 10.2…]₅₃
[8.5…, 9.25…]₅₃
[7.5…, 8.5…]₅₃
[6.5…, 7.5…]₅₃
[5.62…, 6.5…]₅₃
[4.88…, 5.62…]₅₃
[3.0…, 4.75…]₅₃
[2.19…, 3.0…]₅₃
[1.28…, 2.19…]₅₃
[-0.5…, 1.31…]₅₃

julia> ps =[ # from https://discourse.julialang.org/t/root-isolation-of-real-rooted-integer-polynomials/51421/1
                      942438915208811912419937422298363203125
                   164182217245953398816894035758761902846875
                  4584900574568933770264468813466870772155175
                 48995332393110515074735075708247882042540865
                266674183150777010544241114017741621823207005
                852443280934837985352088128423887894438557515
               1738546146302892245736990844587487000484756535
               2381558158813900978436742173305983349418813145
               2262003889258248241081177038309445610985409335
               1516025051068561122302699855213604175083575145
                720810987764796866354279087485114863858085005
                241341213302325116160821160849326697595681275
                 55563538585119205063994483187179167406616375
                  8363912237256094118085070946083688310200625
                   740493466743082745510080711751444519503125
                    29215606371473169285018060091249259296875];

julia> ANewDsc(ps)
There were 15 isolating intervals found:
[-0.01617…, 0.0531…]₂₅₆
[-0.08643…, -0.01617…]₂₅₆
[-0.2285…, -0.08643…]₂₅₆
[-0.3711…, -0.2285…]₂₅₆
[-0.6562…, -0.3672…]₂₅₆
[-0.9531…, -0.6562…]₂₅₆
[-1.234…, -0.9531…]₂₅₆
[-1.84…, -1.25…]₂₅₆
[-2.125…, -1.812…]₂₅₆
[-2.438…, -2.125…]₂₅₆
[-3.0…, -2.44…]₂₅₆
[-3.281…, -3.031…]₂₅₆
[-3.562…, -3.281…]₂₅₆
[-3.875…, -3.562…]₂₅₆
[-4.188…, -3.875…]₂₅₆
```

## Refinement

The `refine_roots` method for the state object uses a bisection method to identify the roots.

## Comparisons

Don't judge the algorithm by its implementation here. This
implementation is not as performant as it could be.

Comparing to some alternatives, we have that the functionality from
Hecke.jl (`Hecke._roots`) is **much** better. (The last example is 33
≈ 0.14s/0.0042 times faster)

However, compared to other alternatives this implementation could be seen as useful:

```
julia> x = variable(Polynomial);

julia> p = -1 + 254*x - 16129*x^2 + x^15;

julia> ANewDsc(coeffs(p))  # ≈ 0.05 seconds;
There were 3 isolating intervals found:
[0.75…, 4.25…]₅₃
[0.00787401574803149653139…, 0.00787401574803149972047…]₂₁₂
[0.0078740157480314937167…, 0.00787401574803149653139…]₂₁₂

julia> filter(isreal, roots(p)) # much faster, but misses two roots with imaginary part ~ 1e-10
1-element Array{Complex{Float64},1}:
 2.1057742291764914 + 0.0im

julia> filter(isreal, AMRVW.roots(Float64.(coeffs(p)))) # using AMRVW. Similarly misses two roots
1-element Array{Complex{Float64},1}:
 2.1057742291764407 + 0.0im

julia> filter(isreal, AMRVW.roots(BigFloat.(coeffs(p)))) # this works here
3-element Vector{Complex{BigFloat}}:
 0.007874015748031494751793842937491860399146218747427882112208862187156603046408902 + 0.0im
 0.007874015748031497374190409031015351762713667398750832139185145345899098824178322 + 0.0im
    2.105774229176482954331883383107195997983004314462374928263620342390986189650864 + 0.0im

julia> filter(isreal, PolynomialRoots.roots(coeffs(p))) # using PolynomialRoots. Misses 2.105?
2-element Array{Complex{Float64},1}:
   0.0078740158234482 + 0.0im
 0.007874015672614794 + 0.0im

julia> IntervalRootFinding.roots(x->p(x), IntervalArithmetic.Interval(0.0, 5.0)) # using IntervalRootFinding, IntervalArithmetic
8-element Array{Root{IntervalArithmetic.Interval{Float64}},1}:
 Root([2.10577, 2.10578], :unique)
 Root([0.00787418, 0.00787422], :unknown)
 Root([0.00787403, 0.00787409], :unknown)
[...]

julia> @syms x::real # using SymPy

julia> @time  rts = sympy.real_roots(p(x)); # correctly identifies 3.
  0.003896 seconds (518 allocations: 13.359 KiB)

julia> sympy.N(rts[2]) # takes a long time! (162s)
0.00787401574803150
```

The algorithm used is a partial implementation of one presented in:

*Computing Real Roots of Real Polynomials ... and now For Real!*
by Alexander Kobel, Fabrice Rouillier, Michael Sagraloff
[arXiv](https://arXiv.org/1605.00410); [DOI:](https://doi.org/10.1145/2930889.2930937).

and

*Computing real roots of real polynomials*
Michael Sagraloff, Kurt Mehlhorn
[DOI:](https://doi.org/10.1016/j.jsc.2015.03.004)

The algorithm relies on Descartes' rule of signs, which gives a bound
on the number of positive real roots of a polynomial, `p(x)`. By
considering the polynomial `p((ax+b)/(x+1))` (a mapping of `[a,b]` to
`[0,∞)` a bound on the number of roots in `[a,b]` can be found. A
simple bisection algorithm for a square-free polynomial would be to
start with an interval large enough to contain all the real roots,
then subdivide at the midpoint throwing out sub intervals which are
found to have no root; saving intervals known to have 1 root, and
repeating the subdivision otherwise. Issues with this are the need for
many subdivisions when clusters of roots are present and the numeric
issues that arise in computing the mapping.

The work of Sagaraloff, Melhorn, Kobel, and Rouillier improves this by
introducing a Newton test and boundary test for rapidly decreasing the
size of an interval when possible, and the ability to use finite
precision arithmetic, instead of exact arithmetic, to compute the
Descartes' bound, in addition to other algorithmic improvements (not
all implemented here).

!!! note
    A square free polynomial can be found through `p/gcd(p, p')`,
    though in practice this calculation is numerically unstable.

!!! note
    This implementation is **much** slower than the `Hecke.roots`
    function provided through `arblib` in `Hecke.jl`, which itself
    says is not competitive with more specialized implementations,
    such as provided by the paper authors
    (http://anewdsc.mpi-inf.mpg.de/). There are several reasons: The
    `mobius_transform!` function is `𝑶(n²)`, and could be `𝑶(n⋅log(n))`
    with more effort; the polynomial evaluation in `admissiblepoint`
    could, similarly, be made more efficient; despite using tricks
    learned from the `MutableArithmetics.jl` package to reduce
    allocations with the `BigFloat` type, there are still *far* too
    many allocations as each time the precision is changed new
    (allocating) values must be created, as the old ones can't be
    modified in place (using `set_prec!` causes segfaults); the
    significant engineering speedups suggested by Kobel, Rouillier,
    and Sagraloff are not implemented; etc.

"""
ANewDsc, refine_roots

function ANewDsc(p::Container{<:Real}; root_bound=root_bound(p), max_depth=32)

    T = BigFloat

    n = length(p) - 1
    n <= 0 && return State{n+1,T,eltype(p)}(𝐈{T}[],  𝐈{T}[], NTuple{n+1,eltype(p)}(p))

    max_depth *= ceil(Int, log2(n+2))

    n′ = sum(divrem(n,2))
    p′ = [i*p[i+1] for i ∈ 1:n]

    m′,M′ = T.(root_bound)
    I = 𝐈(m′, M′)

    Internal = [(I,1,1)]
    Isol = 𝐈{T}[]
    Unresolved = 𝐈{T}[]
    ctr = [0,0,0]
    while !isempty(Internal)

        I, nᴵ, depth = pop!(Internal)
        depth += 1
        L = precision(I)

        # 01-test is decisive
        zerotest(p, I) && continue
        onetest(p, I) && (push!(Isol, I); continue)

        ## pump the brakes, if needed, otherwise can go on and on in a newtonstep
        ## with non-square input
        depth >= 10max_depth && (push!(Unresolved, I); continue)

        # shrink or divide;
        val, J = newtontest(p, p′, I, nᴵ)
        if val
            push!(Internal, (J,nᴵ+1, depth))
            continue
        end

        val, J = boundarytest(p, I, nᴵ)
        if val
            push!(Internal, (J,nᴵ,depth))
            continue
        end

        # bisection
        Iₗ,Iᵣ = bisect_interval(p, I, n′, L)
        nᴵ = max(0, nᴵ-1)
        append!(Internal, ((Iₗ,nᴵ,depth), (Iᵣ,nᴵ,depth)))

    end

    return State(Isol, Unresolved, NTuple{n+1,eltype(p)}(p))

end
