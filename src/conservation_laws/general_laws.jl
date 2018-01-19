
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
