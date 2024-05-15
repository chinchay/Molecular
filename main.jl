include("atom.jl")
include("sample.jl")

using .AtomModule
using .SampleModule



atom1 = Atom(0, 1, 2)
sample = Sample([atom1])

for a in sample.vAtoms
    println("atom1 x = ", a.x, " ", a.y, " ", a.z)
end