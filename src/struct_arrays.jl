using .StructArrays: StructArray, components

function _subtypes(type::Type)
  out = Any[]
  _subtypes!(out, type)
end

function _subtypes!(out, type::Type)
  if !isabstracttype(type)
    push!(out, type)
  else
    foreach(T -> _subtypes!(out, T), subtypes(type))
  end
  out
end


# Use lispy tuple programming to iterate over all fields of a `StructArray`.
# Use ugly `@eval` to resolve ambiguities arising from existing `mul!` definitions.
for Dtype in _subtypes(AbstractDerivativeOperator)
  @eval Base.@propagate_inbounds function mul!(dest::A, D::$Dtype, u::A) where {A<:StructArray}
    iterate_mul!(Tuple(components(dest)), D, Tuple(components(u)))
  end

  @eval Base.@propagate_inbounds function iterate_mul!(dest::NTuple{N}, D::$Dtype, u::NTuple{N}) where N
    mul!(first(dest), D, first(u))
    iterate_mul!(Base.tail(dest), D, Base.tail(u))
  end


  @eval Base.@propagate_inbounds function mul!(dest::A, D::$Dtype, u::A, α) where {A<:StructArray}
    iterate_mul!(Tuple(components(dest)), D, Tuple(components(u)), α)
  end

  @eval Base.@propagate_inbounds function iterate_mul!(dest::NTuple{N}, D::$Dtype, u::NTuple{N}, α) where N
    mul!(first(dest), D, first(u), α)
    iterate_mul!(Base.tail(dest), D, Base.tail(u), α)
  end


  @eval Base.@propagate_inbounds function mul!(dest::A, D::$Dtype, u::A, α, β) where {A<:StructArray}
    iterate_mul!(Tuple(components(dest)), D, Tuple(components(u)), α, β)
  end

  @eval Base.@propagate_inbounds function iterate_mul!(dest::NTuple{N}, D::$Dtype, u::NTuple{N}, α, β) where N
    mul!(first(dest), D, first(u), α, β)
    iterate_mul!(Base.tail(dest), D, Base.tail(u), α, β)
  end
end

iterate_mul!(::Tuple{}, D::AbstractDerivativeOperator, ::Tuple{}, args...) = nothing
