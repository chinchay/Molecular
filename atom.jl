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
    _r :: Position
    _v :: Velocity
    _a :: Acceleration
end

function Atom()
    r = Position(0, 0, 0)
    v = Velocity(1, 0, 0)
    a = Acceleration(0, 0, 0)
    return Atom(r, v, a)
end

function move(atom::Atom)
    atom._r.x += DT * atom._v.x
    atom._r.y += DT * atom._v.y
    atom._r.z += DT * atom._v.z
end

end