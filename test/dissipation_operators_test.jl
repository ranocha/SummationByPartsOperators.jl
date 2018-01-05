using Base.Test
using SummationByPartsOperators

D_test_list = (MattssonSvärdShoeybi2008(), Mattsson2014(), MattssonAlmquistCarpenter2014Extended())
Di_test_list = (MattssonSvärdNordström2004(),)


# Test symmetry and eigenvalues
for source_D in D_test_list, source_Di in Di_test_list, acc_order in 2:2:8, T in (Float32,Float64)
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    der_order = 1

    D = try
        derivative_operator(source_D, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end

    if D != nothing
        @inferred mass_matrix(D)
        H = mass_matrix(D)

        for order in 2:2:8
            Di = try
                dissipation_operator(source_Di, D, order)
            catch
                nothing
            end
            if Di != nothing
                HDi = H*full(Di)
                @test norm(HDi - HDi') < 10*eps(T)
                @test minimum(real, eigvals(HDi)) > -10*eps(T)
            end
        end
    end
end
