import .DiffEqCallbacks: SavingCallback

function SavingCallback(semi::AbstractSemidiscretization; kwargs...)
  T = eltype(semi.derivative)

  save_func = (u,t,integrator) -> integrate(u->ScalarIntegralQuantities(u,u^2),
                                              u, integrator.f.f)
  saved_values = SavedValues(T, ScalarIntegralQuantities{T})
  SavingCallback(save_func, saved_values; kwargs...)
end
