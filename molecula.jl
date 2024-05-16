module MoleculaModule

import ..AtomModule: Atom

export Molecula

mutable struct Molecula
    vAtoms :: Vector{Atom}

end

end