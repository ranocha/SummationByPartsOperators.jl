using Test, SummationByPartsOperators

for T in (Float32, Float64), split_form in (Val{true}(), Val{false}())
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    tspan = (zero(T), one(T))
    
    # Fourier
    let D = fourier_derivative_operator(xmin, xmax, N)
        Di = dissipation_operator(TadmorWaagan2012Standard(), D)
        semidisc = CubicPeriodicSemidiscretisation(D, Di, split_form)
        ode = semidiscretise(u0func, semidisc, tspan)
        du = similar(ode.u0)
        semidisc(du, ode.u0, nothing, first(tspan))
    end

    # Periodic FD
    for acc_order in 2:2:8
        D = periodic_derivative_operator(1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semidisc = CubicPeriodicSemidiscretisation(D, Di, split_form)
        ode = semidiscretise(u0func, semidisc, tspan)
        du = similar(ode.u0)
        semidisc(du, ode.u0, nothing, first(tspan))
    end

    # Legendre
    let D = legendre_derivative_operator(xmin, xmax, N)
        Di = nothing
        semidisc = CubicNonperiodicSemidiscretisation(D, Di, split_form, zero, zero)
        ode = semidiscretise(u0func, semidisc, tspan)
        du = similar(ode.u0)
        semidisc(du, ode.u0, nothing, first(tspan))
    end

    # SBP FD
    for acc_order in 2:2:8
        D = derivative_operator(MattssonSvärdNordström2004(), 1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semidisc = CubicNonperiodicSemidiscretisation(D, Di, split_form, zero, zero)
        ode = semidiscretise(u0func, semidisc, tspan)
        du = similar(ode.u0)
        semidisc(du, ode.u0, nothing, first(tspan))
    end
end
