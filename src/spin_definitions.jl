eye(n) = sparse(1.0I, n, n)

# generic spin_half
sz_half = sparse(Float64[.5 0; 0 -.5])
sp_half = sparse(Float64[ 0 1; 0  0])
sm_half = sparse(Float64[ 0 0; 1  0])
I2 = eye(2)

# U1 symmetric spin_half
sz_half_U1sym = SymTensor(0, (STLeg(+1, [0, 1], [1, 1]), STLeg(-1, [0, 1], [1, 1])),
                          [(0, 0), (1, 1)],
                          [-.5*ones(1,1), .5*ones(1,1)])
