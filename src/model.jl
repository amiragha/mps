struct Point1D
    x :: Float64
end

struct Point2D
    x :: Float64
    y :: Float64
end

struct Bond1D{T<:Number}
    one :: Int
    two :: Int
    offx :: Int
    t   :: T
end

struct Bond2D{T<:Number}
    one :: Int
    two :: Int
    offx :: Int
    offy :: Int
    t   :: T
end

struct UnitCell1D{T<:Number}
    n :: Int
    sites :: Vector{Point1D}
    vs :: Vector{Point1D}
    bonds :: Vector{Bond1D{T}}
    mus :: Vector{T}
end

struct UnitCell2D{T<:Number}
    n :: Int
    sites :: Vector{Point2D}
    vs :: Vector{Point2D}
    bonds :: Vector{Bond2D{T}}
    mus :: Vector{T}
end

function chainunitcell(t::T, mu::T=zero(T)) where {T<:Number}
     uc = UnitCell1D(
        1,
        [Point1D(0,0)],
        [Point1D(+1)],
        [Bond2D(1, 1, +1, t)],
        [mu])
    uc
end

function squareunitcell(tx::T, ty::T, mu::T=zero(T)) where {T<:Number}
    uc = UnitCell2D(
        1,
        [Point2D(0,0)],
        [Point2D(0,1), Point2D(1,0)],
        [Bond2D(1, 1, +1, 0, tx), Bond2D(1, 1, 0, +1, ty)],
        [mu])
    uc
end

function triangular_unitcell(t1::T, t2::T, mu::T=zero(T)) where {T<:Number}
    uc = UnitCell2D(
        1,
        [Point2D(0,0)],
        [Point2D(0.5,sin(pi/3)), Point2D(1,0)],
        [Bond2D(1, 1, +1, 0, t2), Bond2D(1, 1, 0, +1, t1), Bond2D(1, 1, -1, +1, t1)],
        [mu])
    uc
end

##TODO: combine this stuff into a function that works for every dimension.
function makemodel(
    uc::UnitCell2D{T},
    lx::Int,
    ly::Int;
    boundary::Symbol=:OBC) where{T<:Number}

    nsites = uc.n*ly*lx
    H = zeros(T, nsites, nsites)
    for (y,x) in Iterators.product(1:ly,1:lx)
        l = (x-1) * ly * uc.n + (y-1) * uc.n
        for n in 1:uc.n
            H[l+n,l+n] = uc.mus[n]
        end
        for b in uc.bonds
            index1 = l+b.one
            if 0 < y+b.offy <= ly && 0 < x+b.offx <= lx
                index2 = l+b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                #println("$x, $y => ($index1, $index2)")
                H[index1, index2] += -b.t
                H[index2, index1] += -conj(b.t)
            elseif boundary == :PBCX && 0 < y+b.offy <= ly
                x+b.offx > lx ? lnew = l - lx * ly * uc.n : lnew = l + lx * ly * uc.n
                index2 = lnew + b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                #println("wrap x, $y => ($index1, $index2)")
                H[index1, index2] += -b.t
                H[index2, index1] += -conj(b.t)
            elseif !(boundary in [:OBC, :PBCX])
                error("unrecognized boundary condition : ", boundary)
            end
        end
    end
    H
end
