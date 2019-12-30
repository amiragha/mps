abstract type AbstractGroup end
abstract type AbstractAbelianGroup <: AbstractGroup end

abstract type AbstractCharge end
struct Trivial <: AbstractCharge end

"U1Charge for U(1) abelian group"
struct U1Charge <: AbstractCharge
    charge :: Int
end
const U1 = U1Charge

@inline Base.inv(c::U1Charge) = U1Charge(-c.charge)
@inline Base.:+(c1::U1Charge, c2::U1Charge) = U1Charge(c1.charge + c2.charge)
@inline Base.:+(c::U1Charge) = c
@inline Base.:-(c1::U1Charge, c2::U1Charge) = U1Charge(c1.charge - c2.charge)
@inline Base.:-(c::U1Charge) = inv(c)
@inline Base.zero(::Type{U1Charge}) = U1Charge(0)
@inline Base.isless(c1::U1Charge, c2::U1Charge) = isless(c1.charge, c2.charge)
@inline Base.isequal(c1::U1Charge, c2::U1Charge) = isequal(c1.charge, c2.charge)
convert(::Type{U1Charge}, s::Int) = U1Charge(s)

function Base.show(io::IO, c::AbstractCharge)
    show(c.charge)
    nothing
end
