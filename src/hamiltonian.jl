module hamiltonian
    using ITensors
    export TFIM_Hamiltonian, XXZ_Hamiltonian, ground_state, excited_state
    export Heisenberg_Hamiltonian, WZW_2_2_Hamiltonian
    export H_dict
    export dim_dict


    function TFIM_Hamiltonian(N,
                            Jz::Float64=-1.,
                            h::Float64=-1.)
        pbc = true
        sites = siteinds("S=1/2",N)
        ampo = AutoMPO()

        # Input operator terms which define
        # a Hamiltonian matrix, and convert
        # these terms to an MPO tensor network
        for j=1:N-1
          add!(ampo, 4 * Jz, "Sz", j, "Sz", j+1)
        end
        if (pbc)
            add!(ampo, 4 * Jz, "Sz", N, "Sz", 1)
        end
        for j=1:N
          add!(ampo, 2 * h, "Sx", j)
        end

        H = MPO(ampo,sites)
        return H, sites
    end


    function WZW_2_2_Hamiltonian(N)
        sites = siteinds("S=1",N)
        ampo = AutoMPO()
        for j=1:N-1

          # S dot S
          add!(ampo, 2, "S+", j, "S-", j+1)
          add!(ampo, 2, "S-", j, "S+", j+1)
          add!(ampo, 4, "Sz", j, "Sz", j+1)

          # - (S dot S) ^2
          add!(ampo, -1, "S+", j, "S-", j+1, "S+", j, "S-", j+1)
          add!(ampo, -1, "S+", j, "S-", j+1, "S-", j, "S+", j+1)
          add!(ampo, -1, "S-", j, "S+", j+1, "S-", j, "S+", j+1)
          add!(ampo, -1, "S-", j, "S+", j+1, "S+", j, "S-", j+1)

          add!(ampo, -2, "S+", j, "S-", j+1, "Sz", j, "Sz", j+1)
          add!(ampo, -2, "S-", j, "S+", j+1, "Sz", j, "Sz", j+1)
          add!(ampo, -2, "Sz", j, "Sz", j+1, "S+", j, "S-", j+1)
          add!(ampo, -2, "Sz", j, "Sz", j+1, "S-", j, "S+", j+1)

          add!(ampo, -4, "Sz", j, "Sz", j+1, "Sz", j, "Sz", j+1)
        end
        # PBC
        add!(ampo, 2, "S+", N, "S-", 1)
        add!(ampo, 2, "S-", N, "S+", 1)
        add!(ampo, 4, "Sz", N, "Sz", 1)

        add!(ampo, -1, "S+", N, "S-", 1, "S+", N, "S-", 1)
        add!(ampo, -1, "S+", N, "S-", 1, "S-", N, "S+", 1)
        add!(ampo, -1, "S-", N, "S+", 1, "S-", N, "S+", 1)
        add!(ampo, -1, "S-", N, "S+", 1, "S+", N, "S-", 1)

        add!(ampo, -2, "S+", N, "S-", 1, "Sz", N, "Sz", 1)
        add!(ampo, -2, "S-", N, "S+", 1, "Sz", N, "Sz", 1)
        add!(ampo, -2, "Sz", N, "Sz", 1, "S+", N, "S-", 1)
        add!(ampo, -2, "Sz", N, "Sz", 1, "S-", N, "S+", 1)

        add!(ampo, -4, "Sz", N, "Sz", 1, "Sz", N, "Sz", 1)

        # Convert these terms to an MPO tensor network

        H = MPO(ampo,sites)
        return H, sites
    end


    function Heisenberg_Hamiltonian(N)
        # Input operator terms which define a Hamiltonian
        sites = siteinds("S=1",N)
        ampo = AutoMPO()
        for j=1:N-1
          add!(ampo, 1, "S+", j, "S-", j+1)
          add!(ampo, 1, "S-", j, "S+", j+1)
          add!(ampo, 2, "Sz", j, "Sz", j+1)
        end
        # PBC
        add!(ampo, 1, "S+", N, "S-", 1)
        add!(ampo, 1, "S-", N, "S+", 1)
        add!(ampo, 2, "Sz", N, "Sz", 1)

        # Convert these terms to an MPO tensor network

        H = MPO(ampo,sites)
        return H, sites
    end

    function XXZ_Hamiltonian(N,
                             J::Float64=1.,
                             Delta::Float64=.5)
        # Input operator terms which define a Hamiltonian
        sites = siteinds("S=1/2",N)
        ampo = AutoMPO()
        for j=1:N-1
          add!(ampo, 4 * J/2., "S+", j, "S-", j+1)
          add!(ampo, 4 * J/2., "S-", j, "S+", j+1)
          add!(ampo, 4 * Delta, "Sz", j, "Sz", j+1)
        end
        # PBC
        add!(ampo, 4 * J/2., "S+", N, "S-", 1)
        add!(ampo, 4 * J/2., "S-", N, "S+", 1)
        add!(ampo, 4 * Delta, "Sz", N, "Sz", 1)

        # Convert these terms to an MPO tensor network

        H = MPO(ampo,sites)
        return H, sites
    end

    function ground_state(H, sites, sweeps)
        psi0 = randomMPS(sites)

        # Run the DMRG algorithm, returning energy
        # (dominant eigenvalue) and optimized MPS
        energy, psi = dmrg(H,psi0, sweeps,outputlevel=1)
        println("Final energy = $energy")
        energy, psi
    end

    function excited_state(psis, H, sites, sweeps)
        # Takes in previously found wave functions and finds orthogonal excited state
        psi0 = randomMPS(sites)

        energy, psi = dmrg(H,psis, psi0, sweeps,outputlevel=1, weight=100)
        energy, psi
    end

    dim_dict = Dict("TFIM" => 2,       # spin 1/2
                    "XXZ" => 2,        # spin 1/2
                    "Heisenberg" => 3, # spin 1
                    "WZW" => 3         # spin 1
                    )
    H_dict = Dict("TFIM" => TFIM_Hamiltonian,
                  "XXZ" => XXZ_Hamiltonian,
                  "Heisenberg" => Heisenberg_Hamiltonian,
                  "WZW" => WZW_2_2_Hamiltonian
                  )
end
