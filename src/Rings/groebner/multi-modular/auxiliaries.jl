function _generate_primes_set(first::Int=2^28, last::Int=-1)
    return PrimesSet(first, last)
end

function _get_next_primes(ps::PrimesSet, prev_prime::Int64, nr::Int=1)
    p   = prev_prime
    return [p = iterate(ps, p)[1] for i in 1:nr]
end

function _modular_groebner_run(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool, prime::Int)
    Qt    = base_ring(I)
    R     = GF(prime)
    Rt, t = PolynomialRing(R, [string(s) for s = symbols(Qt)], cached = false)
    Ip    = ideal(Rt, [Rt(f) for f in gens(I)])
    
    if ordering == degrevlex(gens(base_ring(I)))
        return f4(Ip, la_option=44, complete_reduction=complete_reduction)
    else
	    return groebner_basis(Ip, ordering = ordering, complete_reduction = complete_reduction)
    end
end

function _lift_modular_basis(Gp::Vector{gfp_mpoly}, ZR::FmpzMPolyRing)
    @debug "Lifting step"
    return map(f->lift(ZR, f), Gp)
end

function _modular_groebner_loop(I::MPolyIdeal, ordering::MonomialOrdering,
		complete_reduction::Bool, primes::Vector{Int64}, ZR::FmpzMPolyRing)

    G = Vector{Vector{gfp_mpoly}}(undef, length(primes))
    H = Vector{Vector{fmpz_mpoly}}(undef, length(primes))
    Threads.@threads for i in 1:length(primes)
        G[i]  = _modular_groebner_run(I, ordering, complete_reduction, primes[i])
        H[i]  = _lift_modular_basis(G[i], ZR)
    end

    return H
end

function _induce_rational_reconstruction(f::fmpz_mpoly, d::fmpz, b::Bool; parent=1)
  g = MPolyBuildCtx(parent)
  for (c, v) in zip(coefficients(f), exponent_vectors(f))
    fl, r, s = Hecke.rational_reconstruction(c, d)
    if !fl
      return false, finish(g)
    end
    push_term!(g, r//s, v)
  end
  return true, finish(g)
end

function _induce_crt(f::fmpz_mpoly, d::fmpz, g::fmpz_mpoly, p::fmpz, b::Bool)
  mu = MPolyBuildCtx(parent(f))
  for i=1:length(f)
    e = exponent_vector(f, i)
    @assert e == exponent_vector(g, i)
    push_term!(mu, crt(coeff(f, i), d, coeff(g, i), p, b), e)
  end
  return finish(mu), d*p
end

function _rational_reconstruction!(G::Vector{Vector{fmpz_mpoly}},
        primes::Vector{Int}, basis::Vector{Vector{fmpq_mpoly}},
        intermediate::Vector{Vector{fmpz_mpoly}}, run::Vector{Int},
        d::Vector{fmpz}, TR::MPolyRing)
    for i in 1:length(primes)
        if basis[1] != []
            h = G[i]
            new_idx = [any(x -> !iszero(GF(primes[i])(x)),
                coefficients(map_coefficients(QQ, h[j], parent = TR) - basis[1][j]))
                for j in 1:length(intermediate[1])]
            fl = !any(new_idx)
            if !fl
                for j in 1:length(intermediate[1])
                    if (new_idx[j])
                        intermediate[1][j], _ =
                        _induce_crt(intermediate[1][j], d[1], h[j], fmpz(primes[i]), true)
                    end
                end
                d[1]  *=  fmpz(primes[i])
                fl    =   true
                for j in 1:length(intermediate[1])
                    if new_idx[j]
                        fl, basis[1][j] = _induce_rational_reconstruction(
                                   	  	intermediate[1][j], d[1], false, parent = TR)
                        fl || break
                    end
                end
                run[1]  = 2
            else
                d[1]    *=  fmpz(primes[i])
                run[1]  -=  1
            end
        else
            d[1]  *=  fmpz(primes[i])
            intermediate[1] = G[i]
            for f = intermediate[1]
                fl, fQ  = _induce_rational_reconstruction(f, d[1], false, parent = TR)
                fl || break
                push!(basis[1], fQ)
            end
            for j in length(basis[1])+1:length(intermediate[1])
                push!(basis[1], TR(0))
            end
        end
    end
end
