using Test, SummationByPartsOperators

@testset "Var. lin. adv." for T in (Float32, Float64), split_form in (Val{true}(), Val{false}())
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    tspan = (zero(T), one(T))
    afunc = one

    # Legendre
    let D = legendre_derivative_operator(xmin, xmax, N)
        Di = nothing
        semi = VariableLinearAdvectionNonperiodicSemidiscretization(D, Di, afunc, split_form, zero, zero)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))
    end

    # SBP FD
    for acc_order in 2:2:8
        D = derivative_operator(MattssonSvärdNordström2004(), 1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semi = VariableLinearAdvectionNonperiodicSemidiscretization(D, Di, afunc, split_form, zero, zero)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))
    end
end
