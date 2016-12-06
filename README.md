Geometric Algebra for Julia
======

This is work-in-progress straight port of [Fontijne's reference implementation][impl] of geometric
algebra utilities to Julia. The code is right now mostly unperformant and non-idiomatic, but I try to
improve it whenever I can.

This project follows the license of the original implementation, GPL2.

[impl]: http://www.geometricalgebra.net/reference_impl.html "Fontijne's implementation in Java"

Things to do
------

A list of things to do, in no particular order:

 - Rewrite in idiomatic Julia
 - Implement general_inverse and utility functions for Multivectors
 - Support named basis vectors
 - Maybe a builtin support for popular geometry models (conformal, hyperbolic, homogeneous etc)
 - A comprehensive test suite
 - cos() and sin()
 - Multivector types
 - meet and join (possibly using the new method described [here][newpaper]
 - A comprehensive benchmark to track performance improvements/regressions
 - Travis-ci support?
 - Write documentation
 - Use BitArray for bitmaps instead of integers?

[newpaper]: http://www.geometricalgebra.net/downloads/fontijne_agacse2008_fact_join_blades.pdf "Fontijne's paper"

Quick start
-----------

julia> include("GeoAlg.jl")

julia> using GeoAlg

julia> e0 = basisvector(0)

1.0

julia> e1 = basisvector(1)

1.0*e1

julia> e2 = basisvector(2)

1.0*e2

julia> e3 = basisvector(3)

1.0*e3

julia> scalarproduct(e1,e1)

1.0

julia> scalarproduct(e1,e2)

0.0

julia> e1 + 2 * e2 + e3

1.0*e1 + 2.0*e2 + 1.0*e3

julia> e1^e2

1.0*e1^e2

julia> e1^e1

0

julia> e1^e2

1.0*e1^e2

julia> (e0 + e1 + e1^e3)^e2

1.0*e2 + 1.0*e1^e2 - 1.0*e1^e2^e3

julia> (e0 + e1 + e1^e3)^e2^e3

1.0*e2^e3 + 1.0*e1^e2^e3
