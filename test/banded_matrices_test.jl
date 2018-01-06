using Base.Test
using SummationByPartsOperators, BandedMatrices

D_test_list = (MattssonSvärdShoeybi2008(), Mattsson2014(), MattssonAlmquistCarpenter2014Extended(),
                MattssonAlmquistCarpenter2014Optimal())
Di_test_list = (MattssonSvärdNordström2004(),)

for T in (Float32, Float64), acc_order in (2,4,6,8), D_source in D_test_list, Di_source in Di_test_list
    xmin = zero(T)
    xmax = 5*one(T)
    N = 101
    D_serial = try
        derivative_operator(D_source, 1, acc_order, xmin, xmax, N, Val{:serial}())
    catch
        nothing
    end
    D_serial == nothing && continue

    D_full   = full(D_serial)
    D_sparse = sparse(D_serial)
    D_banded = BandedMatrix(D_serial)
    x = grid(D_serial)
    u = cospi.(x)
    dest1 = zeros(u)
    dest2 = zeros(u)

    @test BandedMatrices.isbanded(D_serial) == BandedMatrices.isbanded(D_banded)
    @test bandwidth(D_serial, 1) == bandwidth(D_banded, 1)
    @test bandwidth(D_serial, 2) == bandwidth(D_banded, 2)

    A_mul_B!(dest1, D_serial, u)
    A_mul_B!(dest2, D_full, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    A_mul_B!(dest2, D_sparse, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    A_mul_B!(dest2, D_banded, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))

    Di_serial = try
        derivative_operator(Di_source, 1, acc_order, xmin, xmax, N, Val{:serial}())
    catch
        nothing
    end
    Di_serial == nothing && continue

    Di_full   = full(Di_serial)
    Di_sparse = sparse(Di_serial)
    Di_banded = BandedMatrix(Di_serial)

    @test BandedMatrices.isbanded(Di_serial) == BandedMatrices.isbanded(Di_banded)
    @test bandwidth(Di_serial, 1) == bandwidth(Di_banded, 1)
    @test bandwidth(Di_serial, 2) == bandwidth(Di_banded, 2)

    A_mul_B!(dest1, Di_serial, u)
    A_mul_B!(dest2, Di_full, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    A_mul_B!(dest2, Di_sparse, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    A_mul_B!(dest2, Di_banded, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
end
