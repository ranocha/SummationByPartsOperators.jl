module SummationByPartsOperatorsStructArraysExt

if isdefined(Base, :get_extension)
  using StructArrays: StructArray, components
else
  using ..StructArrays: StructArray, components
end

using SummationByPartsOperators: AbstractDerivativeOperator
import SummationByPartsOperators: mul!

using InteractiveUtils: subtypes

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
  @eval Base.@propagate_inbounds function mul!(dest::StructArray, D::$Dtype, u::StructArray)
    iterate_mul!(Tuple(components(dest)), D, Tuple(components(u)))
    return dest
  end

  @eval Base.@propagate_inbounds function mul!(dest::StructArray, D::$Dtype, u::StructArray, α)
    iterate_mul!(Tuple(components(dest)), D, Tuple(components(u)), α)
    return dest
  end

  @eval Base.@propagate_inbounds function mul!(dest::StructArray, D::$Dtype, u::StructArray, α, β)
    iterate_mul!(Tuple(components(dest)), D, Tuple(components(u)), α, β)
    return dest
  end
end

Base.@propagate_inbounds function iterate_mul!(dest::NTuple{N, Any},
                                               D::AbstractDerivativeOperator,
                                               u::NTuple{N, Any}, args...) where N
  mul!(first(dest), D, first(u), args...)
  iterate_mul!(Base.tail(dest), D, Base.tail(u), args...)
end

iterate_mul!(::Tuple{}, D::AbstractDerivativeOperator, ::Tuple{}, args...) = nothing

end # module
