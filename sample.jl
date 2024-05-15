module SampleModule

import ..AtomModule: Atom

export Sample

mutable struct Sample
    vAtoms :: Vector{Atom}

end

end