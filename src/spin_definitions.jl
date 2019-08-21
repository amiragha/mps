eye(n) = sparse(1.0I, n, n)

"""
    spinoperators(S)

For spin `S` returns the operators (`Sz`,`Sp`, `Sm`)
"""

function spinoperators(S)
    d = Int(round(2*S + 1))
    ms = [S-n+1 for n in 1:d]
    amps = [sqrt((S+m)*(S-m+1)) for m in ms[1:d-1]]

    Sz = diagm(0 => ms)
    Sp = diagm(1 => amps)
    Sm = diagm(-1 => amps)

    Sz, Sp, Sm
end

# generic spin_half
const sz_half = sparse(Float64[.5 0; 0 -.5])
const sp_half = sparse(Float64[ 0 1; 0  0])
const sm_half = sparse(Float64[ 0 0; 1  0])
const I2 = eye(2)

# U1 symmetric spin_half
const sz_half_U1sym = SymTensor(0,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(0, 0), (1, 1)],
                                [-.5*ones(1,1), .5*ones(1,1)])

const sp_half_U1sym = SymTensor(+1,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(1, 0)],
                                [ones(1,1)])

const sm_half_U1sym = SymTensor(-1,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(1, 0)],
                                [ones(1,1)])
