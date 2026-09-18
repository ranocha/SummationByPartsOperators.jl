module SummationByPartsOperatorsForwardDiffExt

if isdefined(Base, :get_extension)
    using ForwardDiff: Dual, Partials
else
    using ..ForwardDiff: Dual, Partials
end

using SummationByPartsOperators: SummationByPartsOperators,
                                 FourierDerivativeOperator,
                                 FourierPolynomialDerivativeOperator,
                                 FourierRationalDerivativeOperator,
                                 PeriodicRationalDerivativeOperator
import SummationByPartsOperators: mul!

# `Dual`s and `Partials` are isbits structs storing their components
# contiguously. Thus, vectors of such values can be reinterpreted as matrices
# of the underlying native numbers, which allows LoopVectorization.jl to
# vectorize the kernels of `FastMode` across the components.
# See https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
function SummationByPartsOperators.reinterpreted_components(::Type{D}) where {Tag, V, N,
                                                                              D <:
                                                                              Dual{Tag,
                                                                                   V, N}}
    SummationByPartsOperators.reinterpreted_components(D, Val{N + 1}(), V)
end

function SummationByPartsOperators.reinterpreted_components(::Type{P}) where {N, V,
                                                                              P <:
                                                                              Partials{N,
                                                                                       V}}
    SummationByPartsOperators.reinterpreted_components(P, Val{N}(), V)
end

# FFTW.jl cannot handle `Dual`s and `Partial`s.
# Thus, we need to specialize the behavior here. It would be even better to
# use the same approach for `Dual`s and an arbitrary number of partials, but
# that doesn't work since FFTW.jl cannot handle non-unit strides.
for Dtype in [FourierDerivativeOperator,
    FourierPolynomialDerivativeOperator,
    FourierRationalDerivativeOperator,
    PeriodicRationalDerivativeOperator]
    @eval Base.@propagate_inbounds function mul!(dest::AbstractVector{Partials{1, T}},
                                                 D::$Dtype,
                                                 u::AbstractVector{Partials{1, T}}) where {T}
        _dest = reinterpret(reshape, T, dest)
        _u = reinterpret(reshape, T, u)
        mul!(_dest, D, _u)
        return dest
    end
end

end # module
