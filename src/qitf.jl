"""
    qitf_hamiltonian(lx, g, J)

generate the quantum Ising in a transverse field Haimltonian that is H
= J ∑_i s^z_i S^z_{i+1} + g ∑_i s^x_i where `lx` is the size of the
system, `J` ising exchange amplitude and `g` is the strength of the
transverse field. The phase transition happens at `g=J` which is the
default value.

"""
function qitf_hamiltonian(lx::Int, g::Float64=1.0, J::Float64=-1.0;
                          envelope::Vector{Float64} = ones(lx),
                          boundary::Symbol=:OBC)
    sx = sp_half + sm_half
    sz = 2*sz_half
    bond_term = J * kron(sz, sz) + g/2. * (kron(sx, I2) + kron(I2, sx))
    Hmat = envelope[1] * bond_term
    if lx > 2
        for i=2:lx-1
            Hmat = kron(Hmat, I2) + envelope[i] * kron(eye(2^(i-1)), bond_term)
        end
    end
    Hmat += envelope[1] * (g/2 * kron(sx, eye(2^(lx-1))) + g/2 * kron(eye(2^(lx-1)), sx))
    if boundary == :OBC
        return Hmat
    elseif boundary == :PBC
        Hmat += envelope[lx] * J * kron(kron(sz, eye(2^(lx-2))), sz)
        return Hmat
    else
        error("unrecognized boundary condition :", boundary)
    end
end

function qitf_bondtensor(g::Float64=1.0, J::Float64=-1.0)
    sx = sp_half + sm_half
    sz = 2*sz_half

    bond_term = Matrix(J * kron(sz, sz) + g/2. * (kron(sx, I2) + kron(I2, sx)))
    return reshape(bond_term, 2,2,2,2)
end

function qitf_energy_exact(g::Float64=1.0, J::Float64=-1.0)
        e_exact, e_exact_err = quadgk(k -> -sqrt(1+g^2-2*g*cos(k))/pi, 0, pi)
    e_exact
end
