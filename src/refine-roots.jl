## -----
## refinement

# For [a,b] an isolating interval for a square free polynomial
# p(a) and p(b) have different signs.
# this finds a refined isolating interval [a′, b′] with |b′ - a′| < 2^{-L}
# if found within 64 refinement steps.
function refine_interval(p, a, b, L=min(53,maximum(precision, (a,b))))
    n = length(p) - 1
    n′ = sum(divrem(n,2))

    p′ = [i*p[i+1] for i ∈ 1:n]

    ϵ = epsᴸ(L)
    I = 𝐈(a,b)
    N = 1
    bisect_depth = 0

    while true

        if bisect_depth > 64
            return I
        end
        
        val, J = newtontest(p, p′, I, N)   
        if val
            diff(J) < ϵ  && return J
            λ = diff(I)/diff(J)
            I = J
            N += 1
            continue
        end

        val, J = boundarytest(p, I, N)
        if val
            diff(J) < ϵ  && return J
            I = J
            continue
        end
        

        Iₗ, Iᵣ = bisect_interval(p, I, n′, L)
        
        val = signtest(p, Iₗ)
        if val == nothing
            val = signtest(p, Iᵣ)
            val != nothing && (val = !val)
        end
        
        if val == nothing
            dₗ,dᵣ = diff(Iₗ), diff(Iᵣ)
            if dₗ <= dᵣ
                dₗ ≤ ϵ && return Iₗ
            else
                dᵣ ≤ ϵ && return Iᵣ
            end
            error("I did not cleanly divide")
        end


        I = val ? Iₗ : Iᵣ
        diff(I) ≤ ϵ && return I
        N = max(0, N-1)
        bisect_depth += 1
    end
end

# refine intervals to size `L`, return midpoint
refine_roots(st::State) = midpoint.(refine_interval(st.p, I...) for I ∈ st)



