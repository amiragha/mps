abstract type AbstractGroup end
abstract type AbstractAbelianGroup <: AbstractGroup end

abstract type AbstractCharge end
struct Trivial <: AbstractCharge end
abstract type AbstractNonTrivialCharge end


"U1Charge for U(1) abelian group"
struct U1Charge <: AbstractCharge
    charge :: Int
end
const U1 = U1Charge

@inline Base.inv(c::U1Charge) = U1Charge(-c.charge)
@inline Base.:*(a::Int, c::U1Charge) = U1Charge(a*c.charge)
@inline Base.:+(c1::U1Charge, c2::U1Charge) = U1Charge(c1.charge + c2.charge)
@inline Base.:+(c::U1Charge, a::Int) = U1Charge(c.charge + a)
@inline Base.:+(c::U1Charge) = c
@inline Base.:-(c1::U1Charge, c2::U1Charge) = U1Charge(c1.charge - c2.charge)
@inline Base.:-(c::U1Charge, a::Int) = U1Charge(c.charge - a)
@inline Base.:-(c::U1Charge) = inv(c)
@inline Base.div(c::U1Charge, a::Int) = div(c.charge, a)
@inline Base.zero(::Type{U1Charge}) = U1Charge(0)
@inline Base.isodd(c::U1Charge) = isodd(c.charge)
@inline Base.isless(c1::U1Charge, c2::U1Charge) = isless(c1.charge, c2.charge)
@inline Base.isequal(c1::U1Charge, c2::U1Charge) = isequal(c1.charge, c2.charge)
@inline Base.isequal(c1::U1Charge, c2::Int) = isequal(c1.charge, c2)
@inline Base.length(::U1Charge) = 1
@inline Base.iterate(c::U1Charge) = (c, nothing)
convert(::Type{U1Charge}, s::Int) = U1Charge(s)

function Base.show(io::IO, c::AbstractCharge)
    print(io, c.charge)
    nothing
end
