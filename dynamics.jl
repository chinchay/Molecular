
# include("position.jl")
# include("velocity.jl")
# include("constants.jl")

module DynamicsModule

include("atom.jl")
include("constants.jl")

import .ConstantsModule
import .ConstantsModule: s, ϵ, s6, s12
using LinearAlgebra # functions: norm

import ..AtomModule: Atom

export updateHarmonicAccelaration, updateLennarJonesAcceleration

function computeDistance(atom1::Atom, atom2::Atom)
    dR21 = atom2._r - atom1._r
    distance = norm(dR21)
    return dR21, distance
end

function computeLennardJonesForce(distance::Float64)
    repulsive = s12 / (distance ^ 13)
    attractive = 0.5 * s6 / (distance ^ 7)
    force = 48 * ϵ * (repulsive - attractive)
    return force
end

function updateLennarJonesAcceleration(atom1::Atom, atom2::Atom)
    dR21, distance = computeDistance(atom1, atom2)
    absForce = computeLennardJonesForce(distance)

    u21 = dR21 / distance
    force21 = absForce * u21
    force12 = -force21
    atom2._a = force21 / atom2.mass
    atom1._a = force12 / atom1.mass
end

function updateHarmonicAccelaration(atom1::Atom, atom2::Atom)
    dR21, distance = computeDistance(atom1, atom2)

    factor = -1
    force21 = factor * dR21
    force12 = -force21

    atom2._a = force21 / atom2.mass
    atom1._a = force12 / atom1.mass

end

end