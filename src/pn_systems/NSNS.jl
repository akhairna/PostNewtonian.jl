"""
    NSNS{T, PNOrder}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There are two additional fields `Λ₁`
and `Λ₂` holding the (constant) tidal-coupling parameters of the neutron stars.  See also
[`BHNS`](@ref).
"""
struct NSNS{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST
    Λ₁::NT
    Λ₂::NT

    NSNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function NSNS(;
        M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Λ₁, Λ₂, Φ=0, PNOrder=typemax(Int), kwargs...
    )
        NT, ST, PNOrder, state = prepare_system(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        return new{NT,ST,PNOrder}(state, convert(NT, Λ₁), convert(NT, Λ₂))
    end
    function NSNS(state; Λ₁, Λ₂, PNOrder=typemax(Int))
        @assert length(state) == 14
        ST, PNOrder = typeof(state), prepare_pn_order(PNOrder)
        NT = eltype(state)
        return new{NT,ST,PNOrder}(state, convert(NT, Λ₁), convert(NT, Λ₂))
    end
end
const BNS = NSNS

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `NSNS` systems.
state(pnsystem::NSNS) = pnsystem.state
function symbols(::Type{<:NSNS})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ, :Λ₁, :Λ₂, )
end
function ascii_symbols(::Type{<:NSNS})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi, :Lambda1, :Lambda2,)
end
