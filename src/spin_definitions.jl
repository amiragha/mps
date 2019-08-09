eye(n) = sparse(1.0I, n, n)

# generic spin_half
const sz_half = sparse(Float64[.5 0; 0 -.5])
const sp_half = sparse(Float64[ 0 1; 0  0])
const sm_half = sparse(Float64[ 0 0; 1  0])
const I2 = eye(2)

# U1 symmetric spin_half
const sz_half_U1sym = SymTensor(0, (STLeg(+1, [0, 1], [1, 1]), STLeg(-1, [0, 1], [1, 1])),
                          [(0, 0), (1, 1)],
                          [-.5*ones(1,1), .5*ones(1,1)])
