
"""
    semidiscretize(u0func, semi::AbstractSemidiscretization, tspan)

Apply the semidiscretization `semi` to the initial data given by `u0func`
and return an `ODEProblem` with time span `tspan`.
"""
function semidiscretize(u0func, semi::AbstractSemidiscretization, tspan)
    u0 = compute_coefficients(u0func, semi.derivative)
    ode = ODEProblem(semi, u0, tspan)
end


function evaluate_coefficients(u, semi::AbstractSemidiscretization)
    evaluate_coefficients(u, semi.derivative)
end


function integrate(func, u, semi::AbstractSemidiscretization)
    integrate(func, u, semi.derivative)
end



@auto_hash_equals struct ScalarIntegralQuantities{T} <: FieldVector{2,T}
    mass::T
    energy::T
end



function FilterCallback(filter::AbstractFilter)
    condition = (t, u, integrator) -> true
    affect!(integrator) = filter(integrator.u)

    DiscreteCallback(condition, affect!, save_positions=(false,false))
end
