
# include("position.jl")
# include("velocity.jl")
# include("constants.jl")

module DynamicsModule

using LinearAlgebra # functions: norm

include("atom.jl")

import ..AtomModule: Atom

export updateAccelarationLJ

function updateAccelarationLJ(atom1::Atom, atom2::Atom)
    dR21 = atom2._r - atom1._r
    # dR12 = -dR21

    distance = norm(dR21)
    u21 = dR21 / distance
    # u12 = -u21

    factor = -1
    force21 = factor * dR21
    force12 = -force21

    atom2._a = force21 / atom2.mass
    atom1._a = force12 / atom1.mass

end

end