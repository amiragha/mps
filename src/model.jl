struct Point2D
    x :: Float64
    y :: Float64
end

struct Bond2D{T<:Number}
    one :: Int
    two :: Int
    offx :: Int
    offy :: Int
    t   :: T
end

struct UnitCell2D{T<:Number}
    n :: Int
    sites :: Vector{Point2D}
    vs :: Vector{Point2D}
    bonds :: Vector{Bond2D{T}}
    mus :: Vector{T}
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

function makemodel(uc::UnitCell2D{T}, lx::Int, ly::Int; boundary::Symbol=:OBC) where{T<:Number}
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
            elseif boundary == :PBC
                error("boundary not supported yet!")
            end
        end
    end
    H
end
