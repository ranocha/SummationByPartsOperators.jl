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
    for acc_order = 2:2:8
        D = derivative_operator(MattssonSvärdNordström2004(), 2, acc_order, xmin, xmax, N)
        semi = WaveEquationNonperiodicSemidiscretization(D, left_bc, right_bc)
        for compact in (true, false)
            show(IOContext(devnull, :compact => compact), semi)
        end
        ode = semidiscretize(du0func, u0func, semi, tspan)
        ddu = similar(ode.u0.x[1])
        semi(ddu, ode.u0.x..., nothing, first(tspan))
    end
end
