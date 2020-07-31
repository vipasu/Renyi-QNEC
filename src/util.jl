module util
    using LinearAlgebra
    using ITensors
    using HDF5
    export vN_entropy, get_ket, extract_vectors, relative_entropy, sandwiched_renyi_divergence
    export calculate_relative_entropies, calculate_srds
    export calculate_entropies
    export save_MPS, load_MPS, get_num_wf

    """
    Calculate von Neumann entropy given an array of eigenvalues.
    """
    function vN_entropy(spec)
        entries = [-s * log(s) for s in spec if s > 0]
    #     println(entries)
        return sum(entries)
    end

    """
    Returns part of an MPS up to site "pos".
    """
    function get_ket(M, pos)
        orthogonalize!(M, pos) # equivalent to "position" in cpp code
        ket = M[1]
        for i=2:pos
            ket *= M[i]
        end
        return ket
    end

    """
    Returns the spectral decomposition of reduced density matrix
    """
    function extract_vectors(psi, b::Int64, localdim=2)
        ket = get_ket(psi, b)
        svd_result = svd(ket,tuple([ket.inds[i] for i = 1:b]...,), full=true, cutoff=1e-30)
        S = svd_result.S
        #     println(S)
        n = localdim^b     #     println(n)
        num_elems = size(svd_result.U.store)[1]
        num_columns = num_elems÷n
        U = reshape(svd_result.U.store, (n,num_columns))
        eigen_vals = S.store.^2
        return U, eigen_vals
    end

    """
    Calculate the relative entropy between two states (given spectral decomposition)
    """
    function relative_entropy(V, P, W, Q)
        # V and W are (arrays of) vectors
        # P and Q are their spectra
        # Assumes the dictionary form for now
        # return S(ρ|σ); ρ=VV^†P, σ=WW^†Q

        rho_diag = Float64[]# diagonal rho elements in the basis of sigma
        for w in eachcol(W)
            temp = 0
            for (block, v) in enumerate(eachcol(V))
                temp += sum(v.*w)^2 * P[block]
            end
            push!(rho_diag, temp)
        end
    #     println(rho_diag)
        total =  vN_entropy(rho_diag) - vN_entropy(collect(values(P)))
        # TODO check for support not contained

        for (p, q) in zip(rho_diag, collect(values(Q)))
            if p > 0
                total += p * (log(p) - log(q))
            end
        end
        total
    end

    """
    Calculate the sandwiched renyi divergence between two states
    (given spectral decomposition) and parameter n.
    """
    function sandwiched_renyi_divergence(V, P, W, Q, n)
        num_v = size(V)[2]
        num_w = size(W)[2]
        WdagV = zeros(num_w, num_v)
        for (i, w) in enumerate(eachcol(W))
            for (j, v) in enumerate(eachcol(V))
                WdagV[i, j] = w' * v
            end
        end
        Qn =  Q .^((1-n)/n)
        sandwich_vals = real.(eigen(WdagV * Diagonal(P) * WdagV' * Diagonal(Qn)).values)
        # raise an issue if complex part is large
        # floor the negative eigenvalues if they are below cutoff
            # should this go here or up ahead?
        sandwich_vals = filter(x -> (x > 0) .& (abs(x) > 1e-60), sandwich_vals)
    #     println(real.(sandwich_vals))
        return 1/(n-1) * log(sum(sandwich_vals.^n))
    end

    """
    Calculates the vN entropy of an MPS for many bipartitions.
    """
    function calculate_entropies(psi, N, localdim)
        entropies = Float64[]
        for i=1:2 *N÷3
            V, P = extract_vectors(psi, i, localdim);
            push!(entropies, vN_entropy(P))
        end
        entropies
    end


    """
    Calculates the relative entropy for two states for many bipartitions.
    """
    function calculate_relative_entropies(phi, psi, N, localdim)
        # assumes psi is the vacuum
        s_rel = Float64[]
        for i=1:2 *N÷3
    #         println(i)
            W, Q = extract_vectors(psi, i, localdim);
            V, P = extract_vectors(phi, i, localdim);
            push!(s_rel, relative_entropy(V,P,W,Q))
        end
        s_rel
    end

    """
    Calculates the sandwiched renyi divergence between two states for many bipartitions.
    """
    function calculate_srds(phi, psi, N, n, localdim)
        # assumes psi is the vacuum
        srd = Float64[]
        for i=1:2 *N÷3
    #         println(i)
            W, Q = extract_vectors(psi, i, localdim);
            V, P = extract_vectors(phi, i, localdim);
            push!(srd, sandwiched_renyi_divergence(V,P,W,Q, n))
        end
        srd
    end

    """
    Helper function to load a saved MPS
    """
    function load_MPS(fname, desc)
        fi = h5open("$fname.h5", "r")
        rmps = read(fi, desc, MPS)
        close(fi)
        return rmps
    end

    """
    Helper function to save an MPS to a datafile with key desc.
    """
    function save_MPS(psi, fname, desc)
        fo = h5open("$fname.h5", "cw")
        write(fo, desc, psi)
        close(fo)
    end

    """
    Routine to get the number of wavefunctions saved in a datafile
    """
    function get_num_wf(fname)
        fi = h5open("$fname.h5", "r")
        name_list = names(fi)
        num_wf = length([n for n in name_list if occursin("wf", n)])
        close(fi)
        return num_wf
    end

end
