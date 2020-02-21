"""
    qitf_hamiltonian(lx, g, J)

generate the quantum Ising in a transverse field Haimltonian that is H
= J ∑_i Z_i Z_{i+1} + g ∑_i X_i where `lx` is the size of the
system, `J` ising exchange amplitude and `g` is the strength of the
transverse field. The phase transition happens at `g=J` which is the
default value.

"""
function qitf_hamiltonian(lx::Int, g::Float64=1.0, J::Float64=-1.0;
                          envelope::Vector{Float64} = ones(lx),
                          boundary::Symbol=:OBC)
    sz, sp ,sm = spinoperators(1/2)
    X = sp + sm
    Z = 2 * sz
    bond_term = J * (Z ⊗ Z) + g/2. * (X ⊗ I(2) + I(2) ⊗ X)
    Hmat = envelope[1] * bond_term
    if lx > 2
        for i=2:lx-1
            Hmat = Hmat ⊗ I(2) + envelope[i] * (I(2^(i-1)) ⊗ bond_term)
        end
    end

    # Should the below line be added only for periodic?! Nahhh...
    Hmat += envelope[1] * g/2 * (X ⊗ I(2^(lx-1)) + I(2^(lx-1)) ⊗ X)

    if boundary == :OBC
        return Hmat
    elseif boundary == :PBC
        Hmat += envelope[lx] * J * (Z ⊗ I(2^(lx-2)) ⊗ Z)
        return Hmat
    else
        throw("unrecognized boundary condition : $boundary")
    end
end

function qitf_bondtensor(g::Float64=1.0, J::Float64=-1.0)
    sz, sp ,sm = spinoperators(1/2)
    X = sp + sm
    Z = 2 * sz
    bond_term = J * (Z ⊗ Z) + g/2. * (X ⊗ I(2) + I(2) ⊗ X)
    reshape(Matrix(bond_term), 2, 2, 2, 2)
end

function qitf_energy_exact(g::Float64=1.0, J::Float64=-1.0)
    e_exact, e_exact_err = quadgk(k -> -sqrt(1+g^2-2*g*cos(k))/pi, 0, pi)
    e_exact
end
