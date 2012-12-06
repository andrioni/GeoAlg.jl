Geometric Algebra for Julia
======

This is work-in-progress straight port of [Fontijne's reference implementation][impl] of geometric
algebra utilities to Julia. The code is right now mostly unperformant and non-idiomatic, but I try to
improve it whenever I can.

This project follows the license of the original implementation, GPL2.

[impl]: http://www.geometricalgebra.net/reference_impl.html "Fontijne's implementation in Java"

Things to do
======

A list of things to do, in no particular order:

 - Rewrite as a module/package
 - Support different metrics
 - Maybe a builtin support for popular geometry models (conformal, hyperbolic, homogeneous etc)
 - A comprehensive test suite
 - cos() and sin()
 - Multivector types
 - meet and join (possibly using the new method described [here][newpaper]

[paper]: http://www.geometricalgebra.net/downloads/fontijne_agacse2008_fact_join_blades.pdf "Fontijne's paper"
