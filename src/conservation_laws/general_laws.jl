
"""
    semidiscretise(u0func, semidisc::AbstractSemidiscretisation, tspan)

Apply the semidiscretisation `semidisc` to the 
"""
function semidiscretise(u0func, semidisc::AbstractSemidiscretisation, tspan)
    u0 = compute_coefficients(u0func, semidisc.derivative)
    ode = ODEProblem(semidisc, u0, tspan)
end
