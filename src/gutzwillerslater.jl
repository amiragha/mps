# Gutzwiller projected slater determinant wavefunction
# TODO: Does Gutzwiller need to carry the states?
mutable struct GutzwillerSlater{T}
    states :: Matrix{T}
    configuration :: Vector{Int}

    uplist :: Vector{Int}
    dnlist :: Vector{Int}

    upslater :: DetMatrix
    dnslater :: DetMatrix
end

"""
    random_gutzwiller_half(states)

returns a Gutzwiller wavefunction with randomly chosen half up spins
and half down spins.  The input `states` is a matrix with spatial
orbitals as columns.  The number of orbitals(columns) has to be
exactly half of the number of sites(rows).

"""
function random_gutzwiller_half(states::Matrix)
    num_sites = size(states, 1)
    num_states = size(states, 2)

    # even lattice size and half-filling
    @assert 2*num_states == num_sites

    ### TODO: can these be simplified in one function?
    configuration = random_ones(num_sites, num_states)
    uplist = findall(x->x==1, configuration)
    dnlist = findall(x->x==0, configuration)

    ### TODO: should be possible to slice them from states, directly!
    ## NOTE: while the states(orbitals) are originally given on cols,
    ## here we put states on rows and positions on columns. Notice the
    ## order of the for expression!

    upslater = DetMatrix([states[x,k] for k=1:num_states, x=uplist])
    dnslater = DetMatrix([states[x,k] for k=1:num_states, x=dnlist])

    return GutzwillerSlater(states, configuration,
                            uplist, dnlist,
                            upslater, dnslater)
end

"""
    propose_step!(wf)

propose and, if accepted, update a step for the Gutzwiller
wavefunction. The process is to choose randomly one up and down spin
fermion and swap their positions. The acceptance probability is the
minimum of 1 and square of the determinant ratio.

returns `true` if update the `wf`, and `false` otherwise.
"""
function propose_step!(wf::GutzwillerSlater)
    ### TODO: split these operations into multiple functions.
    # randomly choose an up and down
    ## NOTE: `up/dn_index` are random indeces in `up/dnlist` pointing to
    ## positions, namely, `up/dn_site` in the wf.configuration!
    up_index = rand(1:length(wf.uplist))
    dn_index = rand(1:length(wf.dnlist))
    up_site = wf.uplist[up_index]
    dn_site = wf.dnlist[dn_index]

    # new columns, i.e. just swap up and dn
    new_up_col = wf.dnslater.matrix[:, dn_index]
    new_dn_col = wf.upslater.matrix[:, up_index]

    # find the probability of the change
    up_detrat = det_ratio_1col(wf.upslater, new_up_col, up_index)
    dn_detrat = det_ratio_1col(wf.dnslater, new_dn_col, dn_index)

    # up_detrat = exact_det_ratio_1col(wf.upslater, new_up_col, up_index)
    # dn_detrat = exact_det_ratio_1col(wf.dnslater, new_dn_col, dn_index)

    if rand(Float64) < min(1, norm(up_detrat * dn_detrat)^2)
        # update the Gutzwiller wavefunction
        wf.configuration[up_site] = 0
        wf.configuration[dn_site] = 1

        wf.uplist[up_index], wf.dnlist[dn_index] = dn_site, up_site

        update_detmatrix_1col!(wf.upslater,
                               new_up_col, up_index, up_detrat)
        update_detmatrix_1col!(wf.dnslater,
                               new_dn_col, dn_index, dn_detrat)

        # exact_update_detmatrix_1col!(wf.upslater,
        #                        new_up_col, up_index, up_detrat)
        # exact_update_detmatrix_1col!(wf.dnslater,
        #                        new_dn_col, dn_index, dn_detrat)

        return true
    else
        return false
    end
end

function check_and_update_gutzwiller!(wf::GutzwillerSlater)
   check_and_update_detmatrix!(wf.upslater)
   check_and_update_detmatrix!(wf.dnslater)
end
