# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

module GeoAlg

import 
    Base.*, Base.^, Base.+, Base.-, Base.~,
    Base.==, Base.<=, Base.>=, Base.<, Base.>,
    Base.cmp, Base.copy, Base.exp, Base.isless, Base.reverse, Base.show,
    Base.string

export
    # Types
    BasisBlade,
    Multivector,
    Metric,

    # Functions
    geometricproduct,   # standard or metric-dependent geometric product, also *
    outerproduct,       # standard outer product, also ^
    innerproduct,       # standard or metric-dependent inner product
    scalarproduct,      # scalar product between multivectors
    grade,              # grade of a blade or multivector
    gradeinversion,     # grade inversion
    cliffordconjugate,  # Clifford conjugate of a blade or multivector
    basisvector,        # n-th basis vector as a multivector
    randomvector,       # a random pure vector 
    randomblade,        # a random blade
    randomversor,       # a random versor, can depend on metric
    scalarpart,         # the scalar (grade 0) part of a multivector
    norm_e,             # euclidean norm of a multivector
    norm_e2,            # squared euclidean norm of a multivector
    extractgrade,       # extract the blades of a certain grade
    versorinverse,      # inverse of a versor
    dual,               # dual of a multivector, can depend on metric
    isnull,             # checks if a multivector is null (or close to it)
    isscalar,           # checks if a multivector is a scalar
    # reverse           # reverses a multivector, also ~

    # Constants
    LEFT_CONTRACTION,
    RIGHT_CONTRACTION,
    HESTENES_INNER_PRODUCT,
    MODIFIED_HESTENES_INNER_PRODUCT

#--------------------------------------------------------------
# Constants
#--------------------------------------------------------------

const LEFT_CONTRACTION = 0
const RIGHT_CONTRACTION = 1
const HESTENES_INNER_PRODUCT = 2
const MODIFIED_HESTENES_INNER_PRODUCT = 3

#--------------------------------------------------------------
# Types
#--------------------------------------------------------------

"
Examples:
julia> t = BasisBlade(0)
1.0

julia> t = BasisBlade(1)
1.0*e1

julia> t = BasisBlade(2)
1.0*e2

julia> t = BasisBlade(1,3.14)^BasisBlade(2)
3.14*e1^e2

The behaviour of this type constructor is not obvious to me.  For example,
some values produce bivectors:

julia> BasisBlade(3)
1.0*e1^e2

julia> BasisBlade(5)
1.0*e1^e3
"
type BasisBlade
    bitmap::Int64
    scale::Float64

    BasisBlade{T<:Number}(b::Int, s::T) = new(b, Float64(s))
    BasisBlade(b::Int) = new(b, 1.0)
    BasisBlade(s::Float64) = new(0, s)
    BasisBlade() = new(0, 0.0)
end

"
Construct a multivector from a BasisBlade.

julia> t = BasisBlade(1,3.14)^BasisBlade(2)
3.14*e1^e2

julia> w = BasisBlade(2)
1.0*e2

julia> z = Multivector(w) + Multivector(t)
1.0*e2 + 3.14*e1^e2
"
type Multivector
    blades::Vector{BasisBlade}
    sorted::Bool

    Multivector() = new(BasisBlade[], false)
    Multivector{T<:Number}(s::T) = new([BasisBlade(Float64(s))], false)
    Multivector(B::Vector{BasisBlade}) = new(B, false)
    Multivector(b::BasisBlade) = new([b], false)
end

type Metric
    matrix::Array{Float64,2}
    eigen::Tuple{Array{Float64,1},Array{Float64,2}}
    inveig::Array{Float64,2}
    metric::Vector{Float64}
    isdiag::Bool
    iseuclidean::Bool
    isantieuclidean::Bool

    function Metric(m::Array{Float64,2})
        matrix = copy(m)
        if !issym(matrix)
            error("the metric matrix must be symmetric")
        end

        eigen = eig(matrix)
        inveig = transpose(eigen[2])
        metric = eigen[1]
        isdiago = isdiag(matrix)
        if !isdiago
            iseuclidean = isantieuclidean = false
        else
            iseuclidean = isantieuclidean = true
            for i = 1:size(matrix, 1)
                if matrix[i,i] != 1.0
                    iseuclidean = false
                end
                if matrix[i,i] != -1.0
                    isantieuclidean = false
                end
            end
        end
        new(matrix, eigen, inveig, metric, isdiago, iseuclidean, isantieuclidean)
    end

    function Metric{T<:Number}(m::Matrix{T})
        matrix = Float64(m)
        Metric(matrix)
    end
end

#--------------------------------------------------------------
# Alternate "constructors"
#--------------------------------------------------------------

"
Return a Multivector representation of a vector.

julia> basisvector(4)
1.0*e4

julia> basisvector(4)^basisvector(2)
-1.0*e2^e4

julia> r = basisvector(4)^basisvector(2)*3
-3.0*e2^e4
"
function basisvector(idx::Int)
    return Multivector(BasisBlade(1 << (idx-1)))
end

function randomvector(dim::Int, scale::Float64)
    result = Array(BasisBlade, dim)
    for i = 1:dim
        result[i] = BasisBlade(1 << (i-1), 2 * scale * (rand() - 0.5))
    end
    return Multivector(result)
end

function randomblade(dim::Int, grade::Int, scale::Float64)
    result = Multivector(2 * scale * (rand() - 0.5))
    for i = 1:grade
        result = result ^ randomvector(dim, scale)
    end
    return result
end

function randomversor(dim::Int, grade::Int, scale::Float64)
    result = Multivector(2 * scale * (rand() - 0.5))
    for i = 1:grade
        result = result * randomvector(dim, scale)
    end
    return result
end

function randomversor(dim::Int, grade::Int, scale::Float64, M::Metric)
    result = Multivector(2 * scale * (rand() - 0.5))
    for i = 1:grade
        result = geometricproduct(result, randomvector(dim, scale), M)
    end
    return result
end

function randomversor(dim::Int, grade::Int, scale::Float64, M::Vector{Float64})
    result = Multivector(2 * scale * (rand() - 0.5))
    for i = 1:grade
        result = geometricproduct(result, randomvector(dim, scale), M)
    end
    return result
end

#--------------------------------------------------------------
# Basic "unary" operations/properties
#--------------------------------------------------------------
grade(a::BasisBlade) = bitcount(a.bitmap)
reverse(a::BasisBlade) = BasisBlade(a.bitmap, minusonepow(div(grade(a) * (grade(a) - 1), 2)) * a.scale);
(~)(a::BasisBlade) = reverse(a)
gradeinversion(a::BasisBlade) = BasisBlade(bitmap, minusonepow(grade(a)) * a.scale)
cliffordconjugate(a::BasisBlade) = BasisBlade(bitmap, minusonepow(div(grade(a) * (grade(a) + 1), 2)) * a.scale);

(*)(a::BasisBlade, b::BasisBlade) = gp_op(a, b, false)
(^)(a::BasisBlade, b::BasisBlade) = gp_op(a, b, true)
function geometricproduct(a::BasisBlade, b::BasisBlade, m::Vector{Float64})
    result = a * b
    bitmap = a.bitmap & b.bitmap
    i = 1
    while bitmap != 0
        if (bitmap & 1) != 0
            result.scale *= m[i]
        end
        i += 1
        bitmap = bitmap >> 1
    end
    return result
end

function dual(A::Multivector, dim::Int)
    I = Multivector(BasisBlade((1 << dim) - 1, 1.0))
    return innerproduct(A, versorinverse(I), LEFT_CONTRACTION)
end

function dual(A::Multivector, M::Metric)
    I = Multivector(BasisBlade((1 << size(M.metric, 1)) - 1, 1.0))
    innerproduct(A, versorinverse(I), M, LEFT_CONTRACTION)
end

function dual(A::Multivector, M::Vector{Float64})
    I = Multivector(BasisBlade((1 << size(M, 1)) - 1, 1.0))
    innerproduct(A, versorinverse(I), M, LEFT_CONTRACTION)
end

function isnull(A::Multivector, epsilon::Float64)
    s = norm_e2(A)
    return (s < epsilon * epsilon)
end

function isnull(A::Multivector)
    simplify(A)
    return size(A.blades, 1) == 0
end

function isscalar(A::Multivector)
    if size(A.blades, 1) == 1
        b = A.blades[1]
        return b.bitmap == 0
    end
    return false
end

function scalarpart(A::Multivector)
    s = 0.0
    for i = 1:size(A.blades)[1]
        b = A.blades[i]
        if b.bitmap == 0
            s += b.scale
        end
    end
    return s
end

function norm_e(A::Multivector)
    s = scalarproduct(A, reverse(A));
    if s < 0.0
        return 0.0
    else
        return sqrt(s)
    end
end


function norm_e2(A::Multivector)
    s = scalarproduct(A, reverse(A));
    if s < 0.0
        return 0.0
    else
        return s
    end
end

function reverse(A::Multivector)
    result = Array(BasisBlade, size(A.blades, 1))
    for i = 1:size(A.blades, 1)
        result[i] = reverse(A.blades[i])
    end
    return Multivector(result)
end
(~)(A::Multivector) = reverse(A)

function gradeinversion(A::Multivector)
    result = Array(BasisBlade, size(A.blades, 1))
    for i = 1:size(A.blades, 1)
        result[i] = gradeinversion(A.blades[i])
    end
    return Multivector(result)
end

function cliffordconjugate(A::Multivector)
    result = Array(BasisBlade, size(A.blades, 1))
    for i = 1:size(A.blades, 1)
        result[i] = cliffordconjugate(A.blades[i])
    end
    return Multivector(result)
end

function extractgrade(A::Multivector, G::Vector{Int})
    maxg = 0
    for i = 1:size(G, 1)
        if G[i] > maxg
            maxg = G[i]
        end
    end

    keep = Array(Bool, maxg + 1)
    for i = 1:size(G, 1)
        keep[G[i]+1] = true
    end

    result = BasisBlade[]
    for i = 1:size(A.blades, 1)
        b = A.blades[i]
        g = grade(b)
        if g > maxg
            continue
        elseif keep[g+1]
            push(result, copy(b))
        end
    end
    return Multivector(result)
end

function versorinverse(A::Multivector)
    R = reverse(A)
    s = scalarproduct(A, R)
    if s == 0.0
        error("Non-invertible multivector")
    end
    return R * (1.0 / s)
end

function versorinverse(A::Multivector, M::Metric)
    R = reverse(A)
    s = scalarproduct(A, R, M)
    if s == 0.0
        error("Non-invertible multivector")
    end
    return R * (1.0 / s)
end

function versorinverse(A::Multivector, M::Vector{Float64})
    R = reverse(A)
    s = scalarproduct(A, R, M)
    if s == 0.0
        error("Non-invertible multivector")
    end
    return R * (1.0 / s)
end

#--------------------------------------------------------------
# Basic binary operations
#--------------------------------------------------------------

function scalarproduct(A::Multivector, B::Multivector)
    return scalarpart(innerproduct(A, B, LEFT_CONTRACTION))
end

function scalarproduct(A::Multivector, B::Multivector, M::Metric)
    scalarpart(innerproduct(A, B, M, LEFT_CONTRACTION))
end

function scalarproduct(A::Multivector, B::Float64, M::Metric)
    scalarpart(innerproduct(A, B, M, LEFT_CONTRACTION))
end

function (+)(A::Multivector, b::Float64)
    result = copy(A.blades)
    push(result, BasisBlade(b))
    return Multivector(simplify(copy(result)))
end
(+)(b::Float64, A::Multivector) = A + b
(+)(b::Int, A::Multivector) = A + Float64(b)
(+)(A::Multivector, b::Int) = A + Float64(b)

function (+)(A::Multivector, B::Multivector)
    result = vcat(A.blades, B.blades)
    return Multivector(simplify(copy(result)))
end

(-)(A::Multivector, b::Float64) = A + (-b)
(-)(b::Float64, A::Multivector) = b + (-1.0 * A)
(-)(A::Multivector, b::Float64) = A + (-Float64(b))
(-)(b::Float64, A::Multivector) = Float64(b) + (-1.0 * A)

function (-)(A::Multivector, B::Multivector)
    result = vcat(A.blades, (-1.0 * B).blades)
    return Multivector(simplify(copy(result)))
end

function geometricproduct(a::BasisBlade, b::BasisBlade, M::Metric)
    A = toeigenbasis(M, a)
    B = toeigenbasis(M, b)
    result = BasisBlade[]
    for i = 1:size(A,1), j = 1:size(B,1)
        push(result, geometricproduct(A[i], B[j], M.metric))
    end
    return tometricbasis(M, simplifybasis(result))
end

function innerproduct(a::BasisBlade, b::BasisBlade, typ::Int)
    innerproductfilter(grade(a), grade(b), a * b, typ)
end

function innerproduct(a::BasisBlade, b::BasisBlade, m::Vector{Float64}, typ::Int)
    innerproductfilter(grade(a), grade(b), geometricproduct(a, b, m), typ)
end

function innerproduct(a::BasisBlade, b::BasisBlade, m::Metric, typ::Int)
    innerproductfilter(grade(a), grade(b), geometricproduct(a, b, M), typ)
end

function (*)(A::Multivector, a::Float64)
    if a == 0.0
        return Multivector()
    end
    result = Array(BasisBlade, size(A.blades)[1])
    for i = 1:size(A.blades)[1]
        b = A.blades[i]
        result[i] = BasisBlade(b.bitmap, b.scale * a)
    end
    return Multivector(result)
end

(*)(a::Float64, A::Multivector) = A * a
(*)(a::Int, A::Multivector) = A * Float64(a)
(*)(A::Multivector, a::Int) = A * Float64(a)

function (*)(A::Multivector, B::Multivector)
    result = Array(BasisBlade, size(A.blades)[1] * size(B.blades)[1])
    k = 1
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result[k] = B1 * B2
            k += 1
        end
    end
    return Multivector(simplify(result))
end

geometricproduct(A::Multivector, a::Float64) = A * a
geometricproduct(a::Float64, A::Multivector) = A * a
geometricproduct(A::Multivector, a::Int) = A * Float64(a)
geometricproduct(a::Int, A::Multivector) = A * Float64(a)
geometricproduct(A::Multivector, B::Multivector) = A * B

function geometricproduct(A::Multivector, B::Multivector, M::Metric)
    result = BasisBlade[]
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result = vcat(result, geometricproduct(B1, B2, M))
        end
    end
    return Multivector(simplify(result))
end

function geometricproduct(A::Multivector, B::Multivector, M::Vector{Float64})
    result = Array(BasisBlade, size(A.blades,1) * size(B.blades,1))
    k = 1
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result[k] = geometricproduct(B1, B2, M)
            k += 1
        end
    end
    return Multivector(simplify(result))
end

outerproduct(A::Multivector, B::Multivector) = A ^ B

function (^)(A::Multivector, B::Multivector)
    result = Array(BasisBlade, size(A.blades)[1] * size(B.blades)[1])
    k = 1
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result[k] = B1 ^ B2
            k += 1
        end
    end
    return Multivector(simplify(result))
end

function innerproduct(A::Multivector, B::Multivector, typ::Int)
    result = Array(BasisBlade, size(A.blades)[1] * size(B.blades)[1])
    k = 1
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result[k] = innerproduct(B1, B2, typ)
            k += 1
        end
    end
    return Multivector(simplify(result))
end

function innerproduct(A::Multivector, B::Multivector, M::Metric, typ::Int)
    result = Array(BasisBlade, size(A.blades)[1] * size(B.blades)[1])
    k = 1
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result[k] = innerproduct(B1, B2, M, typ)
            k += 1
        end
    end
    return Multivector(simplify(result))
end

function innerproduct(A::Multivector, B::Multivector, M::Vector{Float64}, typ::Int)
    result = Array(BasisBlade, size(A.blades)[1] * size(B.blades)[1])
    k = 1
    for i = 1:size(A.blades)[1]
        B1 = A.blades[i]
        for j = 1:size(B.blades)[1]
            B2 = B.blades[j]
            result[k] = innerproduct(B1, B2, M, typ)
            k += 1
        end
    end
    return Multivector(simplify(result))
end

#--------------------------------------------------------------
# Miscellania
#--------------------------------------------------------------

exp(A::Multivector) = exp(A, 12)
function exp(A::Multivector, order::Int)
    A2 = compress(A * A)
    if isnull(A2, 1e-8)
        return A + 1
    elseif isscalar(A2)
        a2 = scalarpart(A2)
        if a2 < 0
            alpha = sqrt(-a2)
            return A * (sin(alpha)/alpha) + cos(alpha)
        else alpha = sqrt(a2)
            return A * (sinh(alpha) / alpha) + cosh(alpha)
        end
    else
        return expseries(A, order)
    end
end

exp(A::Multivector, M::Metric) = exp(A, M, 12)
function exp(A::Multivector, M::Metric, order::Int)
    A2 = compress(geometricproduct(A, A, M))
    if isnull(A2, 1e-8)
        return A + 1
    elseif isscalar(A2)
        a2 = scalarpart(A2)
        if a2 < 0
            alpha = sqrt(-a2)
            return A * (sin(alpha)/alpha) + cos(alpha)
        else alpha = sqrt(a2)
            return A * (sinh(alpha) / alpha) + cosh(alpha)
        end
    else
        return expseries(A, M, order)
    end
end

exp(A::Multivector, M::Vector{Float64}) = exp(A, M, 12)
function exp(A::Multivector, M::Vector{Float64}, order::Int)
    A2 = compress(geometricproduct(A, A, M))
    if isnull(A2, 1e-8)
        return A + 1
    elseif isscalar(A2)
        a2 = scalarpart(A2)
        if a2 < 0
            alpha = sqrt(-a2)
            return A * (sin(alpha)/alpha) + cos(alpha)
        else alpha = sqrt(a2)
            return A * (sinh(alpha) / alpha) + cosh(alpha)
        end
    else
        return expseries(A, M, order)
    end
end

#--------------------------------------------------------------
# Internal use functions
#--------------------------------------------------------------

function expseries(A::Multivector, order::Int)
    scale = 1
    maxi = norm_e(A)
    if maxi > 1.0
        scale <<= 1
    end
    while maxi > 1.0
        maxi = maxi / 2
        scale <<= 1
    end

    scaled = A * (1 / scale)

    result = Multivector(1.0)
    tmp = Multivector(1.0)
    for i = 1:order
        tmp = tmp * scaled * (1/i)
        result = result + tmp
    end

    while scale > 1
        result = result * result
        scale >>>= 1
    end

    return result
end

function expseries(A::Multivector, M::Metric, order::Int)
    scale = 1
    maxi = norm_e(A)
    if maxi > 1.0
        scale <<= 1
    end
    while maxi > 1.0
        maxi = maxi / 2
        scale <<= 1
    end

    scaled = A * (1 / scale)

    result = Multivector(1.0)
    tmp = Multivector(1.0)
    for i = 1:order
        tmp = geometricproduct(tmp, scaled * (1/i), M)
        result = result + tmp
    end

    while scale > 1
        result = geometricproduct(result, result, M)
        scale >>>= 1
    end

    return result
end

function expseries(A::Multivector, M::Vector{Float64}, order::Int)
    scale = 1
    maxi = norm_e(A)
    if maxi > 1.0
        scale <<= 1
    end
    while maxi > 1.0
        maxi = maxi / 2
        scale <<= 1
    end

    scaled = A * (1 / scale)

    result = Multivector(1.0)
    tmp = Multivector(1.0)
    for i = 1:order
        tmp = geometricproduct(tmp, scaled * (1/i), M)
        result = result + tmp
    end

    while scale > 1
        result = geometricproduct(result, result, M)
        scale >>>= 1
    end

    return result
end

function innerproductfilter(ga::Int, gb::Int, R::Vector{BasisBlade}, typ::Int)
    result = []
    for i = 1:size(R)[1]
        B = innerproductfilter(ga, gb, R[i], typ)
        if B.scale != 0.0
            push(result, B)
        end
    end
    return result
end

function innerproductfilter(ga::Int, gb::Int, r::BasisBlade, typ::Int)
    if typ == 0
        if (ga > gb) || (grade(r) != (gb - ga))
            return BasisBlade()
        else
            return r
        end
    elseif typ == 1
        if (ga < gb) || (grade(r) != (ga - gb))
            return BasisBlade()
        else
            return r
        end
    elseif typ == 2
        if (ga == 0) || (gb == 0)
            return BasisBlade()
        end
        # drop through to MODIFIED_HESTENES_INNER_PRODUCT
        if abs(ga - gb) == grade(r)
            return r
        else
            return BasisBlade()
        end
    elseif typ == 3
        if abs(ga - gb) == grade(r)
            return r
        else
            return BasisBlade()
        end
    end
end

function isdiag(A::Matrix)
    if istriu(A) && istril(A)
        return true
    end
    return false
end

function transform(a::BasisBlade, M::Array{Float64,2})
    A = BasisBlade[]
    push(A, BasisBlade(a.scale))

    i = 1
    b = a.bitmap
    while b != 0
        if b & 1 != 0
            tmp = BasisBlade[]
            for j= 1:size(M, 1)
                if M[j,i] != 0
                    m = M[j,i]
                    for k = 1:size(A, 1)
                        push(tmp, A[k] ^ BasisBlade(1<<(j-1), m))
                    end
                end
            end
            A = tmp
        end
        b >>= 1
        i += 1
    end
    return A
end

function toeigenbasis(M::Metric, a::BasisBlade)
    transform(a, M.inveig)
end

function tometricbasis(M::Metric, a::BasisBlade)
    transform(a, M.eigen[2])
end

function tometricbasis(M::Metric, a::Vector{BasisBlade})
    result = BasisBlade[]
    for i = 1:size(a,1)
        tmp = tometricbasis(M, a[i])
        result = vcat(result, tmp)
    end
    return simplifybasis(result)
end

function simplifybasis(A::Vector{BasisBlade})
    if size(A,1) == 0
        return A
    end

    sort!((x,y)->basiscmp(x,y) <= 0, A)
    result = BasisBlade[]
    current = copy(A[1])
    for i = 2:size(A,1)
        b = A[i]
        if b.bitmap == current.bitmap
            current.scale += b.scale
        else
            if current.scale != 0.0
                push(result, current)
            end
            current = copy(b)
        end
    end
    if current.scale != 0.0
        push(result, current)
    end
    return result
end

function canonical_reordering_sign(a::Int, b::Int)
    a >>= 1
    s = 0
    while a != 0
        s += bitcount(a & b)
        a >>= 1
    end
    return ((s & 1) == 0) ? 1.0 : -1.0
end

function bitcount(i::Int)
    i = i - ((i >> 1) & 0x55555555)
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333)
    i = (i + (i >> 4)) & 0x0F0F0F0F
    i = i + (i >> 8)
    i = i + (i >> 16)
    return i & 0x0000003F
end

function gp_op(a::BasisBlade, b::BasisBlade, outer::Bool)
    if outer && ((a.bitmap & b.bitmap) != 0)
        return BasisBlade()
    end

    bitmap = a.bitmap $ b.bitmap
    sig = canonical_reordering_sign(a.bitmap, b.bitmap)

    return BasisBlade(bitmap, sig * a.scale * b.scale)
end

function basiscmp(a::BasisBlade, b::BasisBlade)
    return a.bitmap < b.bitmap ? -1 :
    ((a.bitmap > b.bitmap) ? 1 : cmp(a.scale, b.scale))
end

minusonepow(i::Int) = ((i & 1) == 0) ? 1 : -1

(==)(a::BasisBlade, b::BasisBlade) = (a.bitmap == b.bitmap) && (a.scale == b.scale)

function string(a::BasisBlade)
    result = ""
    i = 1
    b = a.bitmap
    while b != 0
        if (b & 1) != 0
            if length(result) > 0
                result *= "^"
            end
            result *= "e$i"
        end
        b >>= 1
        i += 1
    end

    return length(result) == 0 ? string(a.scale) : string(a.scale) * "*" * result
end

function show(io::IO, a::BasisBlade)
    print(io, string(a))
end

# Multivectors

function (==)(A::Multivector, B::Multivector)
    zer = A - B
    return size(zer.blades,1) == 0
end

function string(A::Multivector)
    if size(A.blades, 1) == 0
        return "0"
    end
    result = ""
    for i = 1:size(A.blades, 1)
        b = A.blades[i]
        S = string(b)
        if i == 1
            result *= S
        elseif S[1] == '-'
            result *= " - " * S[2:end]
        else
            result *= " + " * S
        end
    end
    return result
end

function show(io::IO, a::Multivector)
    print(io, string(a))
end

function copy(a::BasisBlade)
    return BasisBlade(a.bitmap, a.scale)
end

function cmp(a::BasisBlade, b::BasisBlade)
    cmp(a.bitmap, b.bitmap)
end

==(a::BasisBlade, b::BasisBlade) = cmp(a,b) == 0
<=(a::BasisBlade, b::BasisBlade) = cmp(a,b) <= 0
>=(a::BasisBlade, b::BasisBlade) = cmp(a,b) >= 0
<(a::BasisBlade, b::BasisBlade) = cmp(a,b) < 0
>(a::BasisBlade, b::BasisBlade) = cmp(a,b) > 0

isless(a::BasisBlade, b::BasisBlade) = a < b

function simplify(L::Vector{BasisBlade})
    sort!(L)
    prev = 0
    removenull = false
    i = 1
    while i <= size(L, 1)
        curblade = L[i]
        if curblade.scale == 0.0
            del(L, i)
            prev = 0
        elseif (prev != 0) && (prev.bitmap == curblade.bitmap)
            prev.scale += curblade.scale
            del(L, i)
        else
            if (prev != 0) && (prev.scale == 0.0)
                removenull = true
            end
            prev = curblade
            i += 1
        end
    end

    if removenull
        i = 1
        while i <= size(L, 1)
            curblade = L[i]
            if curblade.scale == 0.0
                del(L, i)
            else
                i += 1
            end
        end
    end

    return L
end

function simplify(A::Multivector)
    simplify(A.blades)
    return A
end

function compress(A::Multivector, epsilon::Float64)
    simplify(A)
    maxmag = 0.0
    for i = 1:size(A.blades, 1)
        b = A.blades[i]
        maxmag = max(abs(b.scale), maxmag)
    end
    if maxmag == 0.0
        A.blades = BasisBlade[]
    else
        maxmag = epsilon
        i = 1
        while i <= size(A.blades, 1)
            b = A.blades[i]
            if abs(b.scale) < maxmag
                del(A.blades, i)
            else
                i += 1
            end
        end
    end
    return A
end

compress(A::Multivector) = compress(A, 1e-13)

end # module
