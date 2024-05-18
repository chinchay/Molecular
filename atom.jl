module AtomModule

include("position.jl")
include("velocity.jl")
include("acceleration.jl")
include("constants.jl")

import .PositionModule: Position
import .VelocityModule: Velocity
import .AccelerationModule: Acceleration
import .ConstantsModule
import .ConstantsModule: DT

export Atom

mutable struct Atom
    _r :: Vector # Position
    _v :: Vector # Velocity
    _a :: Vector # Acceleration
    mass :: Float16
end

function Atom()
    r = [0, 0, 0] #Position(0, 0, 0)
    v = [0, 0, 0] #Velocity(1, 0, 0)
    a = [0, 0, 0] #Acceleration(0, 0, 0)
    mass = 1
    return Atom(r, v, a, mass)
end

function Atom(pos::Vector, vel::Vector, acc::Vector, mass::Float16)
    # r = pos # Position(pos[1], pos[2], pos[3])
    # v = vel # Velocity(vel[1], vel[2], vel[3])
    # a = acc # Acceleration(acc[1], acc[2], acc[3])
    return Atom(pos, vel, acc, mass)
end

function move(atom::Atom)
    # dt2 = DT * DT

    # atom._r.x += (DT * atom._v.x) #+ (0.5 * atom._a.x * dt2)
    # atom._r.y += (DT * atom._v.y) #+ (0.5 * atom._a.y * dt2)
    # atom._r.z += (DT * atom._v.z) #+ (0.5 * atom._a.z * dt2)

    # atom._v.x += atom._a.x * DT
    # atom._v.y += atom._a.y * DT
    # atom._v.z += atom._a.z * DT


    atom._r += atom._v * DT
    atom._v += atom._a * DT

end

end