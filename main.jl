include("atom.jl")

using .AtomModule
import .AtomModule: move
using Plots

pos = [0, 0, 0]
vel = [1, 0, 0]
atom1 = Atom(pos, vel)

pos = [4, 0, 0]
vel = [0, 1, 0]
atom2 = Atom(pos, vel)

n = 100
anim = @animate for i = 1:n
    move(atom1)
    move(atom2)
    r1 = atom1._r
    r2 = atom2._r
    scatter([r1.x], [r1.y], [r1.z], legend=false, camera=(30,30), xlim=(-2, 10), ylim=(-2, 10), zlim=(-2, 10))
    scatter!([r2.x], [r2.y], [r2.z], legend=false, camera=(30,30), xlim=(-2, 10), ylim=(-2, 10), zlim=(-2, 10))
end

gif(anim, "molecula.gif", fps=30)
