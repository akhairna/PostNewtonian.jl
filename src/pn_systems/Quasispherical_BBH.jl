"""
    QuasisphericalBBH{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a quasispherical binary black hole system.

The `state` vector here holds the fundamental variables `M₁`, `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`,
with the spins unpacked into three components each, and `R` unpacked into four — for a total
of 13 elements.

Optionally, `Φ` may also be tracked as the 14th element of the `state` vector.  This is just
the integral of the orbital angular frequency `Ω`, and holds little interest for general
systems beyond a convenient description of how "far" the system has evolved.
"""
struct QuasisphericalBBH{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST

    QuasisphericalBBH{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function QuasisphericalBBH(; PNOrder=typemax(Int), kwargs...)
        (NT, ST, PNOrder, state) = prepare_system(QuasisphericalBBH; PNOrder, kwargs...)
        return new{NT,ST,PNOrder}(state)
    end
end
const QuasisphericalBHBH = QuasisphericalBBH

function pack_state(::Type{<:QuasisphericalBBH}; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=0)
    [M₁; M₂; vec(QuatVec(χ⃗₁)); vec(QuatVec(χ⃗₂)); components(Rotor(R)); v; Φ]
end

state(pnsystem::QuasisphericalBBH) = pnsystem.state

function symbols(::Type{<:QuasisphericalBBH})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)
end

function ascii_symbols(::Type{<:QuasisphericalBBH})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
end

# function StaticArrays.SVector(pnsystem::QuasisphericalBBH)
#     return SVector{14,eltype(pnsystem)}(
#         pnsystem.state[1],
#         pnsystem.state[2],
#         pnsystem.state[3],
#         pnsystem.state[4],
#         pnsystem.state[5],
#         pnsystem.state[6],
#         pnsystem.state[7],
#         pnsystem.state[8],
#         pnsystem.state[9],
#         pnsystem.state[10],
#         pnsystem.state[11],
#         pnsystem.state[12],
#         pnsystem.state[13],
#         pnsystem.state[14],
#     )
# end

# TODO the @eval's moved to fundamental_variables.jl

# @testitem "PNSystem constructors" begin
#     using Quaternionic

#     pnA = BBH(;
#         M₁=1.0f0, M₂=2.0f0, χ⃗₁=Float32[3.0, 4.0, 5.0], χ⃗₂=Float32[6.0, 7.0, 8.0], v=0.23f0
#     )
#     @test pnA.state ==
#         Float32[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 1.0; 0.0; 0.0; 0.0; 0.23; 0.0]

#     pnB = BBH(;
#         M₁=1.0f0,
#         M₂=2.0f0,
#         χ⃗₁=Float32[3.0, 4.0, 5.0],
#         χ⃗₂=Float32[6.0, 7.0, 8.0],
#         v=0.23f0,
#         Φ=9.0f0,
#     )
#     @test pnB.state ==
#         Float32[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 1.0; 0.0; 0.0; 0.0; 0.23; 9.0]

#     R = randn(RotorF32)
#     pn1 = BBH(;
#         M₁=1.0f0,
#         M₂=2.0f0,
#         χ⃗₁=Float32[3.0, 4.0, 5.0],
#         χ⃗₂=Float32[6.0, 7.0, 8.0],
#         R=R,
#         v=0.23f0,
#     )
#     @test pn1.state ≈ [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 0.0]

#     pn2 = BBH(;
#         M₁=1.0f0,
#         M₂=2.0f0,
#         χ⃗₁=Float32[3.0, 4.0, 5.0],
#         χ⃗₂=Float32[6.0, 7.0, 8.0],
#         R=R,
#         v=0.23f0,
#         Φ=9.0f0,
#     )
#     @test pn2.state ≈ [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 9.0]

#     pn1.state[end] = 9.0f0
#     @test pn1.state == pn2.state
# end
