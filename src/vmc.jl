mutable struct VMCgutzwiller
    model :: UnitCellQModel

    wavefunction :: GutzwillerSlater

    total_steps :: Int

    n_proposed :: Int
    n_accepted :: Int

    data :: Vector{Float64}
end

function report(sim::VMCgutzwiller)
    println("Finished with ", sim.total_steps)
    println("Number of proposed moves = ", sim.n_proposed)
    println("Number of accepted moves = ", sim.n_accepted)
    n_sites = numofsites(sim.model)
    Dict(zip([(i,j) for i=1:n_sites for j=i+1:n_sites],
             sim.data./sim.n_accepted))
end

function runVMC(model::UnitCellQModel,
                total_steps::Int)
    n_sites = numofsites(model)

    hmat = generatebdg(model)
    symmetrize = eltype(hmat) <: Complex ? Hermitian : Symmetric
    fact = eigen(symmetrize(hmat))

    n_occupied = div(n_sites, 2)
    println(fact.values[n_occupied-2:n_occupied+2])
    e_fermi = fact.values[n_occupied]
    n=1
    while fact.values[n_occupied+n] - e_fermi < 1.e-10
        n=n+1
    end
    states = fact.vectors[:, 1:n_occupied]

    wavefunction = random_gutzwiller_half(states)

    sim = VMCgutzwiller(model, wavefunction,
                        total_steps, 0, 0,
                        zeros(Float64, binomial(n_sites, 2)))

    while sim.n_accepted < sim.total_steps
        sim.n_proposed +=1
        if propose_step!(sim.wavefunction)
            sim.n_accepted += 1

            corrdata = measure(sim.wavefunction, :ZZ)
            sim.data += corrdata

            # display progress every 1000 accepted steps
            if sim.n_accepted % 1000 == 0
                check_and_update_gutzwiller!(sim.wavefunction)
                print("*")
            end
        end
    end
    println()
    report(sim)
end

function measure(wf::GutzwillerSlater, m::Symbol)
    n_sites = size(wf.states, 1)
    n_operators = binomial(n_sites, 2)
    correlations = zeros(Float64, n_operators)
    index = 1
    for i=1:n_sites, j=i+1:n_sites
        ni = wf.configuration[i] - 1/2
        nj = wf.configuration[j] - 1/2
        correlations[index] = ni * nj
        index += 1
    end
    correlations
end
