include("atom.jl")
include("dynamics.jl")

using .AtomModule
import .AtomModule: move
import .DynamicsModule: updateAccelarationLJ
using Plots

pos = [0, 1, 0]
vel = [1, 0, 0]
acc = [0, 0, 0]
mass = 100
atom1 = Atom(pos, vel, acc, mass)

pos = [1, 0, 0]
vel = [0, 1, 0]
acc = [0, 0, 0]
mass = 1
atom2 = Atom(pos, vel, acc, mass)

n = 100
anim = @animate for i = 1:n
    move(atom1)
    move(atom2)
    r1 = atom1._r
    r2 = atom2._r
    scatter( [r1[1]], [r1[2]], [r1[3]], legend=false, camera=(0, 90), xlim=(-1, 10), ylim=(-2, 10), zlim=(-1, 1), xlabel="X", ylabel="Y")
    scatter!([r2[1]], [r2[2]], [r2[3]], legend=false, camera=(0, 90), xlim=(-1, 10), ylim=(-2, 10), zlim=(-1, 1), xlabel="X", ylabel="Y")
    # println("i = ", i, " ", r1[1], " - ", r1[2])

    updateAccelarationLJ(atom1, atom2)
end

gif(anim, "molecula.gif", fps=30)
