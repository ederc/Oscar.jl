function modular_groebner_basis_over_qq(I::MPolyIdeal; ordering::MonomialOrdering=degrevlex(gens(base_ring(I))),  complete_reduction::Bool = false)
    if !haskey(I.gb, ordering)
        _compute_modular_groebner_basis_over_qq(I, ordering, complete_reduction)
    end

    return collect(I.gb[ordering])
end

function _compute_modular_groebner_basis_over_qq(I::MPolyIdeal{fmpq_mpoly}, ordering::MonomialOrdering=degrevlex(gens(base_ring(I))), complete_reduction::Bool=false, verify::Bool=false)
    # generate ring for lifting coefficients
    ZR  = PolynomialRing(ZZ, [string(s) for s = symbols(base_ring(I))], cached = false)[1]

    basis         = Vector{Vector{fmpq_mpoly}}(undef, 1)
    intermediate  = Vector{Vector{fmpz_mpoly}}(undef, 1)

    basis[1] = []
    intermediate[1] = []
    d   = fmpz[1]
    run = Int[2]

    nr_primes = 0

    # initialize primes
    primes_set  = _generate_primes_set()
    last_prime  = iterate(primes_set)[1]

    # run modular computations and liftings until the result stabilizes
    while run[1] > 0
        # get next bunch of primes
        primes      = _get_next_primes(primes_set, last_prime, Threads.nthreads())
        #= primes      = get_next_primes(primes_set, last_prime, Threads.nthreads()) =#
        last_prime  = primes[end]

        # run modular Groebner basis computations and lift results
        G = _modular_groebner_loop(I, ordering, complete_reduction, primes, ZR)
    
        _rational_reconstruction!(G, primes, basis, intermediate, run, d, base_ring(I))
        nr_primes += Threads.nthreads()
    end

    if verify
	    # TODO
    end

    I.gb[ordering]  = BiPolyArray(basis[1], keep_ordering = false, isGB = true)
end

