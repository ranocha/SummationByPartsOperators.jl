using Test
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

D_test_list = (MattssonNordström2004(), MattssonSvärdNordström2004(),
                MattssonSvärdShoeybi2008(), Mattsson2014(),
                MattssonAlmquistCarpenter2014Extended(),
                MattssonAlmquistCarpenter2014Optimal())
Di_test_list = (MattssonSvärdNordström2004(),)
D2var_test_list = (Mattsson2012(),)

for T in (Float32, Float64), acc_order in (2,4,6,8), diss_order in (2,4,6,8), D_source in D_test_list, Di_source in Di_test_list
    xmin = zero(T)
    xmax = 5*one(T)
    N = 101
    D_serial = try
        derivative_operator(D_source, 1, acc_order, xmin, xmax, N, Val{:serial}())
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    D_serial == nothing && continue

    D_full   = Matrix(D_serial)
    D_sparse = sparse(D_serial)
    D_banded = BandedMatrix(D_serial)
    x = grid(D_serial)
    u = cospi.(x)
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    @test BandedMatrices.isbanded(D_serial) == BandedMatrices.isbanded(D_banded)
    @test bandwidth(D_serial, 1) == bandwidth(D_banded, 1)
    @test bandwidth(D_serial, 2) == bandwidth(D_banded, 2)

    mul!(dest1, D_serial, u)
    mul!(dest2, D_full, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=5000*eps(T)), eachindex(u))
    mul!(dest2, D_sparse, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=5000*eps(T)), eachindex(u))
    mul!(dest2, D_banded, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=5000*eps(T)), eachindex(u))


    Di_serial = try
        dissipation_operator(Di_source, D_serial, order=diss_order)
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    Di_serial == nothing && continue

    Di_full   = Matrix(Di_serial)
    Di_sparse = sparse(Di_serial)
    Di_banded = BandedMatrix(Di_serial)

    @test BandedMatrices.isbanded(Di_serial) == BandedMatrices.isbanded(Di_banded)
    @test bandwidth(Di_serial, 1) == bandwidth(Di_banded, 1)
    @test bandwidth(Di_serial, 2) == bandwidth(Di_banded, 2)

    mul!(dest1, Di_serial, u)
    mul!(dest2, Di_full, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    mul!(dest2, Di_sparse, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
    mul!(dest2, Di_banded, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=500*eps(T)), eachindex(u))
end


for T in (Float32, Float64), acc_order in (2,4,6), D2var_source in D2var_test_list
    xmin = zero(T)
    xmax = 5*one(T)
    N = 101
    D2var_serial = try
        var_coef_derivative_operator(D2var_source, 2, acc_order, xmin, xmax, N, one, Val{:serial}())
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    D2var_serial == nothing && continue

    D2var_full   = Matrix(D2var_serial)
    D2var_sparse = sparse(D2var_serial)
    D2var_banded = BandedMatrix(D2var_serial)
    x = grid(D2var_serial)
    u = cospi.(x)
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    @test BandedMatrices.isbanded(D2var_serial) == BandedMatrices.isbanded(D2var_banded)
    @test bandwidth(D2var_serial, 1) == bandwidth(D2var_banded, 1)
    @test bandwidth(D2var_serial, 2) == bandwidth(D2var_banded, 2)

    mul!(dest1, D2var_serial, u)
    mul!(dest2, D2var_full, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=5000*eps(T)), eachindex(u))
    mul!(dest2, D2var_sparse, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=5000*eps(T)), eachindex(u))
    mul!(dest2, D2var_banded, u)
    @test all(i->isapprox(dest1[i], dest2[i], atol=5000*eps(T)), eachindex(u))
end
