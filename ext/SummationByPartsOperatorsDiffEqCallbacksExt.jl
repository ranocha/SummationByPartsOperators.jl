module SummationByPartsOperatorsDiffEqCallbacksExt

if isdefined(Base, :get_extension)
    import DiffEqCallbacks: SavingCallback, SavedValues
else
    import ..DiffEqCallbacks: SavingCallback, SavedValues
end

using SummationByPartsOperators:
                                 SciMLBase, AbstractSemidiscretization,
                                 ScalarIntegralQuantities, integrate

function SavingCallback(semi::AbstractSemidiscretization; kwargs...)
    T = eltype(semi.derivative)

    save_func = (u, t, integrator) -> integrate(u -> ScalarIntegralQuantities(u, u^2),
                                                u,
                                                SciMLBase.unwrapped_f(integrator.f.f))
    saved_values = SavedValues(T, ScalarIntegralQuantities{T})
    SavingCallback(save_func, saved_values; kwargs...)
end

end # module
