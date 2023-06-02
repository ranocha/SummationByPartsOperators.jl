module ForwardDiffExt

if isdefined(Base, :get_extension)
  using ForwardDiff: Partials
else
  using ..ForwardDiff: Partials
end

using SummationByPartsOperators: FourierDerivativeOperator,
                                 FourierPolynomialDerivativeOperator,
                                 FourierRationalDerivativeOperator,
                                 PeriodicRationalDerivativeOperator
import SummationByPartsOperators: mul!

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
