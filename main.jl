include("atom.jl")
include("molecula.jl")

using .AtomModule
using .MoleculaModule



atom1 = Atom(0, 1, 2)
atom2 = Atom(1, 3, 3)
molecula = Molecula([atom1, atom2])


for a in molecula.vAtoms
    println("atom x = ", a.x, " ", a.y, " ", a.z)
end
