
"""
    semidiscretise(u0func, semidisc::AbstractSemidiscretisation, tspan)

Apply the semidiscretisation `semidisc` to the initial data given by `u0func`
and return an `ODEProblem` with time span `tspan`.
"""
function semidiscretise(u0func, semidisc::AbstractSemidiscretisation, tspan)
    u0 = compute_coefficients(u0func, semidisc.derivative)
    ode = ODEProblem(semidisc, u0, tspan)
end


function evaluate_coefficients(u, semidisc::AbstractSemidiscretisation)
    evaluate_coefficients(u, semidisc.derivative)
end


function integrate(func, u, semidisc::AbstractSemidiscretisation)
    integrate(func, u, semidisc.derivative)
end



struct ScalarIntegralQuantities{T} <: FieldVector{2,T}
    mass::T 
    energy::T
end

function DiffEqCallbacks.SavingCallback(semidisc::AbstractSemidiscretisation; kwargs...)
    T = eltype(semidisc.derivative)

    save_func = (u,t,integrator) -> integrate(u->ScalarIntegralQuantities(u,u^2),
                                                u, integrator.f)
    saved_values = SavedValues(T, ScalarIntegralQuantities{T})
    SavingCallback(save_func, saved_values; kwargs...)
end



function FilterCallback(filter::AbstractFilter)
    condition = (t, u, integrator) -> true
    affect!(integrator) = filter(integrator.u)

    DiscreteCallback(condition, affect!, save_positions=(false,false))
end
