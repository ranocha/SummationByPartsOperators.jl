using Base.Test
using SummationByPartsOperators, BandedMatrices

test_list = (MattssonSvärdShoeybi2008(), Mattsson2014(), MattssonAlmquistCarpenter2014Extended(),
    MattssonAlmquistCarpenter2014Optimal())

for T in (Float32, Float64), acc_order in (2,4,6,8), source in test_list
    xmin = zero(T)
    xmax = 5*one(T)
    N = 101
    D_serial = try
        derivative_operator(source, 1, acc_order, xmin, xmax, N, Val{:serial}())
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

    A_mul_B!(dest1, D_serial, u)
    A_mul_B!(dest2, D_full, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    A_mul_B!(dest2, D_sparse, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    A_mul_B!(dest2, D_banded, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
end
