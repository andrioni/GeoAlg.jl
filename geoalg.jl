import Base.*
import Base.reverse
import Base.^
import Base.+
import Base.-
import Base.exp
import Base.cmp
import Base.==
import Base.<=, Base.>=, Base.<, Base.>
import Base.isless
import Base.show
import Base.string

const LEFT_CONTRACTION = 0
const RIGHT_CONTRACTION = 1
const HESTENES_INNER_PRODUCT = 2
const MODIFIED_HESTENES_INNER_PRODUCT = 3

type BasisBlade
    bitmap::Int64
    scale::Float64
end

function BasisBlade(b::Int)
    BasisBlade(b, 1.0)
end

function BasisBlade(s::Float64)
    BasisBlade(0, s)
end

function BasisBlade()
    BasisBlade(0, 0.0)
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

(*)(a::BasisBlade, b::BasisBlade) = gp_op(a, b, false)
(^)(a::BasisBlade, b::BasisBlade) = gp_op(a, b, true)

grade(a::BasisBlade) = bitcount(a.bitmap)
minusonepow(i::Int) = ((i & 1) == 0) ? 1 : -1
reverse(a::BasisBlade) = BasisBlade(a.bitmap, minusonepow(div(grade(a) * (grade(a) - 1), 2)) * a.scale);
gradeinversion(a::BasisBlade) = BasisBlade(bitmap, minusonepow(grade(a)) * a.scale)
cliffordconjugate(a::BasisBlade) = BasisBlade(bitmap, minusonepow(div(grade(a) * (grade(a) + 1), 2)) * a.scale);

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

function innerproductfilter(ga::Int, gb::Int, R, typ::Int)
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

function innerproduct(a::BasisBlade, b::BasisBlade, typ::Int)
    innerproductfilter(grade(a), grade(b), a * b, typ)
end

(==)(a::BasisBlade, b::BasisBlade) = (a.bitmap == b.bitmap) && (a.scale == b.scale)

function string(a::BasisBlade)
    result = ""
    i = 1
    b = a.bitmap
    while b != 0
        if (b & 1) != 0
            if size(result)[1] > 0
                result *= "^"
            end
            result *= "e$i"
        end
        b >>= 1
        i += 1
    end

    return size(result)[1] == 0 ? string(a.scale) : string(a.scale) * "*" * result
end

function show(io, a::BasisBlade)
    print(io, string(a))
end


type Multivector
    blades::Vector{BasisBlade}
    sorted::Bool
end

function Multivector()
    blades = BasisBlade[]
    Multivector(blades, false)
end

function Multivector(B::Vector{BasisBlade})
    blades = B
    Multivector(blades, false)
end

function Multivector(s::Float64)
    Multivector([BasisBlade(s)], false)
end

function Multivector(b::BasisBlade)
    Multivector([b], false)
end

function (==)(A::Multivector, B::Multivector)
    zer = A - B
    return size(zer.blades)[1] == 0
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

function show(io, a::Multivector)
    print(io, string(a))
end

function basisvector(idx::Int)
    return Multivector(BasisBlade(1 << idx))
end

function randomvector(dim::Int, scale::Float64)
    result = Array(BasisBlade, dim)
    for i = 1:dim
        result[i] = BasisBlade(1 << i, 2 * scale * (rand - 0.5))
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

# Depende da métrica
#function randomversor(dim::Int, grade::Int, scale::Float64)
    #randomversor(dim, grande, scale, )

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
(*)(a::Int, A::Multivector) = A * float64(a)
(*)(A::Multivector, a::Int) = A * float64(a)

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

function scalarproduct(A::Multivector, B::Multivector)
    return scalarpart(innerproduct(A, B, LEFT_CONTRACTION))
end

function (+)(A::Multivector, b::Float64)
    #result = Array(size(A.blades)[1] + 1)
    result = copy(A.blades)
    push(result, BasisBlade(b))
    return Multivector(simplify(copy(result)))
end
(+)(b::Float64, A::Multivector) = A + b
(+)(b::Int, A::Multivector) = A + float64(b)
(+)(A::Multivector, b::Int) = A + float64(b)

function (+)(A::Multivector, B::Multivector)
    result = vcat(A.blades, B.blades)
    return Multivector(simplify(copy(result)))
end

(-)(A::Multivector, b::Float64) = A + (-b)
(-)(b::Float64, A::Multivector) = b + (-1.0 * A)
(-)(A::Multivector, b::Float64) = A + (-float64(b))
(-)(b::Float64, A::Multivector) = float64(b) + (-1.0 * A)

function (-)(A::Multivector, B::Multivector)
    result = vcat(A.blades, (-1.0 * B).blades)
    return Multivector(simplify(copy(result)))
end

function compress(A::Multivector, epsilon::Float64)
    simplify(A)
    maxmag = 0.0
    for i = 1:size(A.blades, 1)
        b = blades[i]
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


exp(A::Multivector) = exp(A, 12)
# Falta caso com métrica
function exp(A::Multivector, order::Int)
    A2 = compress(A * A)
    if isnull(A2, 1e-8)
        return A + 1
    elseif isscalar(A2)
        a2 = scalarpart(A2)
        if a2 < 0
            alpha = sqrt(-a2)
            return A * sin(alpha)/alpha + cos(alpha)
        else alpha = sqrt(a2)
            return A * sinh(alpha) / alpha + cosh(alpha)
        end
    else
        return expseries(A, order)
    end
end

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

