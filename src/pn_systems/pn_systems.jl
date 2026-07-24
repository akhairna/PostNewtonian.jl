"""
    PNSystem{NT, ST, PNOrder}

Abstract type for all post-Newtonian systems, including quasicircular and
eccentric binary systems such as `BBH`, `BHNS`, and `NSNS`.

Subtypes of `PNSystem` encode the physical properties of a binary system and its
current state. They serve as inputs to the various [fundamental](@ref
Fundamental-variables) and [derived variables](@ref Derived-variables),
as well as PN expressions and [dynamics](@ref Dynamics) functions.

Each subtype should define a `state` vector containing the fundamental variables
for that system. The parameter `NT` specifies the numeric type of the system,
`ST` is the type of the state vector (for example, `Vector{Float64}`), and
`PNOrder` specifies the order to which PN expansions are carried.

Subtypes should provide appropriate `symbols` and `symbol_index` methods to
allow generic access to their fundamental variables.
"""
abstract type PNSystem{NT,ST<:DenseVector{NT},PNOrder} <: DenseVector{NT} end

"""
    state(pnsystem::PNSystem)

Return the state vector of `pnsystem`, which is a vector of fundamental
variables for the given PN system.

Note that the built-in `PNSystem` subtypes typically store the state directly as
a field and this function returns that field. User-defined subtypes may
implement `state` differently, but should return the corresponding state vector.
"""
function state(::T) where {T<:PNSystem}
    error("`state` is not yet defined for PNSystem subtype `$T`.")
end

Base.vec(pnsystem::PNSystem) = state(pnsystem)

const VecOrPNSystem = Union{AbstractVector,PNSystem}

Base.eltype(::Type{PNT}) where {PNT<:PNSystem} = NT
Base.one(::Type{PNT}) where {PNT<:PNSystem} = one(eltype(PNT))
Base.one(x::T) where {T<:PNSystem} = one(T)
Base.zero(::Type{PNT}) where {PNT<:PNSystem} = zero(eltype(PNT))
Base.zero(x::T) where {T<:PNSystem} = zero(T)
Base.float(::Type{PNT}) where {PNT<:PNSystem} = float(eltype(PNT))
Base.float(x::T) where {T<:PNSystem} = float(T)

### Interfaces: https://docs.julialang.org/en/v1/manual/interfaces
# Iteration
Base.iterate(pnsystem::PNSystem) = iterate(state(pnsystem))
Base.iterate(pnsystem::PNSystem, istate) = iterate(state(pnsystem), istate)
Base.IteratorSize(::Type{T}) where {T<:PNSystem} = Base.HasShape{1}()
Base.length(pnsystem::PNSystem) = length(state(pnsystem))
Base.ndims(pnsystem::PNSystem) = ndims(state(pnsystem))
Base.size(pnsystem::PNSystem) = size(state(pnsystem))
Base.size(pnsystem::PNSystem, dim) = size(state(pnsystem), dim)
Base.IteratorEltype(::Type{T}) where {T<:PNSystem} = Base.HasEltype()
Base.eltype(::Type{<:PNSystem{NT}}) where {NT} = NT
Base.isdone(pnsystem::PNSystem) = Base.isdone(state(pnsystem))
Base.isdone(pnsystem::PNSystem, iterstate) = Base.isdone(state(pnsystem), iterstate)

# Indexing
Base.getindex(pnsystem::PNSystem, i::Int) = getindex(state(pnsystem), i)
Base.setindex!(pn::PNSystem, v, i::Int) = setindex!(state(pn), v, i)
Base.firstindex(pnsystem::PNSystem) = firstindex(state(pnsystem))
Base.lastindex(pnsystem::PNSystem) = lastindex(state(pnsystem))
Base.eachindex(pnsystem::PNSystem) = eachindex(state(pnsystem))

# Abstract arrays
Base.IndexStyle(::Type{T}) where {T<:PNSystem} = Base.IndexLinear()
Base.similar(pnsystem::PNSystem) = similar(state(pnsystem))
Base.axes(pnsystem::PNSystem) = axes(state(pnsystem))

# Strided Arrays
Base.strides(pnsystem::PNSystem) = strides(state(pnsystem))
function Base.unsafe_convert(::Type{Ptr{T}}, A::PNSystem) where {T}
    Base.unsafe_convert(Ptr{T}, state(A))
end
Base.elsize(::Type{<:PNSystem{T}}) where {T} = sizeof(T)
Base.stride(pnsystem::PNSystem, k::Int) = stride(state(pnsystem), k)

"""
    pn_order(pnsystem::PNSystem)

Return the PN order of the given `pnsystem`.

This is a `Rational{Int}` that indicates the order to which the PN expansions
should be carried out.
"""
pn_order(::PNSystem{NT,ST,PNOrder}) where {NT,ST,PNOrder} = PNOrder

"""
    order_index(pnsystem::PNSystem)

Return the order index of the given `pnsystem`.

This is defined as the (one-based) index into an iterable of PN terms starting at 0pN, then
0.5pN, etc.  Specifically, this is defined as `1 + Int(2pn_order(pnsystem))`.
"""
order_index(pn::PNSystem) = 1 + Int(2pn_order(pn))

"""
    max_pn_order

The maximum PN order that can be used without overflowing the `Int` type.
"""
const max_pn_order = (typemax(Int) - 2) // 2

"""
    causes_domain_error!(u̇, p)

Ensure that these parameters correspond to a physically valid set of PN parameters.

If the parameters are not valid, this function should modify `u̇` to indicate that the
current step is invalid.  This is done by filling `u̇` with `NaN`s, which will be detected
by the ODE solver and cause it to try a different (smaller) step size.

Currently, the only check that is done is to test that these parameters result in a PN
parameter v>0.  In the future, this function may be expanded to include other checks.
"""
function causes_domain_error!(u̇, p::PNSystem{NT}) where {NT}
    if p.state[symbol_index(typeof(p), Val(:v))] ≤ 0  # If this is expanded, document the change in the docstring.
        u̇ .= convert(NT, NaN)
        true
    else
        false
    end
end

"""
    prepare_system
"""
function prepare_system(T::Type{<:PNSystem}; PNOrder=typemax(Int), kwargs...)
    state = pack_state(T; kwargs...)
    ST = typeof(state)
    NT = eltype(ST)
    PNOrder = prepare_pn_order(PNOrder)
    return (NT, ST, PNOrder, state)
end

"""
    prepare_pn_order(PNOrder)

Convert the input to a half-integer of type `Rational{Int}`.

If `PNOrder` is larger than `max_pn_order`, it is set to `max_pn_order`, to avoid overflow
when computing the order index.
"""
function prepare_pn_order(PNOrder)
    if PNOrder < max_pn_order
        round(Int, 2PNOrder)//2
    else
        max_pn_order
    end
end

"""
    symbols(pnsystem::PNSystem)
    symbols(::Type{<:PNSystem})
    ascii_symbols(pnsystem::PNSystem)
    ascii_symbols(::Type{<:PNSystem})

Return a Tuple of symbols corresponding to the variables tracked by `pnsystem`, in the order
in which they are stored in the `state` vector.

The `ascii_symbols` function returns those symbols in ASCII form, enabling interaction with
external systems (e.g., Python) that do not support many Unicode symbols.

```jldoctest
julia> using PostNewtonian: BBH

julia> pnsystem = BBH(randn(14); PNOrder=7//2);

julia> symbols(pnsystem)
(:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)

julia> ascii_symbols(pnsystem)
(:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
```
"""
symbols(pnsystem::PNSystem) = symbols(typeof(pnsystem))
function symbols(::Type{T}) where {T<:PNSystem}
    error("`symbols` is not yet defined for PNSystem subtype `$T`.")
end

ascii_symbols(pnsystem::PNSystem) = ascii_symbols(typeof(pnsystem))
function ascii_symbols(::Type{T}) where {T<:PNSystem}
    error("`ascii_symbols` is not yet defined for PNSystem subtype `$T`.")
end

symbol_index(pnsystem::PNSystem, s::Symbol) = symbol_index(typeof(pnsystem), Val(s))
function symbol_index(::Type{T}, ::Val{S}) where {T<:PNSystem,S}
    index = findfirst(y -> y == S, symbols(T))
    if isnothing(index)
        index = findfirst(y -> y == S, ascii_symbols(T))
    end
    if isnothing(index)
        error(
            "Type `$(T)` has no symbol `:$(S)`.\n" *
            "This type's symbols are `$(symbols(T))`.\n" *
            "The ASCII equivalents are `$(ascii_symbols(T))`.\n",
        )
    else
        index
    end
end
