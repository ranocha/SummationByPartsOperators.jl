using Test, SummationByPartsOperators

BCs = (Val(:HomogeneousDirichlet), Val(:HomogeneousNeumann), Val(:NonReflecting))

for T in (Float32, Float64), left_bc in BCs, right_bc in BCs
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    du0func = zero
    tspan = (zero(T), one(T))

    # SBP FD
    for acc_order in 2:2:8
        D = derivative_operator(MattssonSvärdNordström2004(), 2, acc_order, xmin, xmax, N)
        semidisc = WaveEquationNonperiodicSemidiscretisation(D, left_bc, right_bc)
        ode = semidiscretise(du0func, u0func, semidisc, tspan)
        ddu = similar(ode.u0[1])
        semidisc(ddu, ode.u0..., nothing, first(tspan))
    end
end
