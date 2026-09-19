using Test
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

D_test_list = (MattssonNordström2004(), MattssonSvärdNordström2004(),
               MattssonSvärdShoeybi2008(), Mattsson2014(),
               MattssonAlmquistCarpenter2014Extended(),
               MattssonAlmquistCarpenter2014Optimal())
Di_test_list = (MattssonSvärdNordström2004(),)
D2var_test_list = (Mattsson2012(),)

for T in (Float32, Float64), acc_order in (2, 4, 6, 8), diss_order in (2, 4, 6, 8),
    D_source in D_test_list, Di_source in Di_test_list

    xmin = zero(T)
    xmax = 5 * one(T)
    N = 101
    D_serial = try
        derivative_operator(D_source, 1, acc_order, xmin, xmax, N, Val{:serial}())
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    D_serial === nothing && continue

    D_full = @inferred Matrix(D_serial)
    D_sparse = @inferred sparse(D_serial)
    D_banded = @inferred BandedMatrix(D_serial)
    D_banded_new = zero(D_banded)
    @inferred copyto!(D_banded_new, D_serial)
    @test D_banded_new == D_banded

    x = grid(D_serial)
    u = cospi.(x)
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    @test BandedMatrices.isbanded(D_serial) == BandedMatrices.isbanded(D_banded)
    @test bandwidth(D_serial, 1) == bandwidth(D_banded, 1)
    @test bandwidth(D_serial, 2) == bandwidth(D_banded, 2)

    # All matrix representations are assembled from the same `mul!` calls and
    # must thus contain exactly the same coefficients. We deliberately use `==`
    # instead of `isapprox` here since the latter compares the difference to the
    # norm of the whole matrix, which is too coarse to detect dropped
    # coefficients of the boundary closure.
    @test D_full == D_sparse
    @test D_full == D_banded

    mul!(dest1, D_serial, u)
    mul!(dest2, D_full, u)
    @test all(i -> isapprox(dest1[i], dest2[i], atol = 5000 * eps(T)), eachindex(u))
    mul!(dest2, D_sparse, u)
    @test dest1 ≈ dest2
    mul!(dest2, D_banded, u)
    @test dest1 ≈ dest2

    Di_serial = try
        dissipation_operator(Di_source, D_serial, order = diss_order)
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    Di_serial === nothing && continue

    Di_full = @inferred Matrix(Di_serial)
    Di_sparse = @inferred sparse(Di_serial)
    Di_banded = @inferred BandedMatrix(Di_serial)
    Di_banded_new = zero(Di_banded)
    @inferred copyto!(Di_banded_new, Di_serial)
    @test Di_banded_new == Di_banded

    @test BandedMatrices.isbanded(Di_serial) == BandedMatrices.isbanded(Di_banded)
    @test bandwidth(Di_serial, 1) == bandwidth(Di_banded, 1)
    @test bandwidth(Di_serial, 2) == bandwidth(Di_banded, 2)

    # All matrix representations must contain exactly the same coefficients
    @test Di_full == Di_sparse
    @test Di_full == Di_banded

    mul!(dest1, Di_serial, u)
    mul!(dest2, Di_full, u)
    @test all(i -> isapprox(dest1[i], dest2[i], atol = 500 * eps(T)), eachindex(u))
    mul!(dest2, Di_sparse, u)
    @test all(i -> isapprox(dest1[i], dest2[i], atol = 500 * eps(T)), eachindex(u))
    mul!(dest2, Di_banded, u)
    @test all(i -> isapprox(dest1[i], dest2[i], atol = 500 * eps(T)), eachindex(u))
end

# The coefficient function `one` results in the constant coefficient case,
# where the boundary closures can be narrower than in the general case.
# Thus, we need to test variable coefficients as well - otherwise, wrong
# bandwidths remain unnoticed, see
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/344
b_test_list = (one, x -> 1 + x^2 / 10)

for T in (Float32, Float64), acc_order in (2, 4, 6), D2var_source in D2var_test_list,
    b_func in b_test_list

    xmin = zero(T)
    xmax = 5 * one(T)
    N = 101
    D2var_serial = try
        var_coef_derivative_operator(D2var_source, 2, acc_order, xmin, xmax, N, b_func,
                                     Val{:serial}())
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    D2var_serial === nothing && continue

    D2var_full = @inferred Matrix(D2var_serial)
    D2var_sparse = @inferred sparse(D2var_serial)
    D2var_banded = @inferred BandedMatrix(D2var_serial)
    D2var_banded_new = zero(D2var_banded)
    @inferred copyto!(D2var_banded_new, D2var_serial)
    @test D2var_banded_new == D2var_banded

    x = grid(D2var_serial)
    u = cospi.(x)
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    @test BandedMatrices.isbanded(D2var_serial) == BandedMatrices.isbanded(D2var_banded)
    @test bandwidth(D2var_serial, 1) == bandwidth(D2var_banded, 1)
    @test bandwidth(D2var_serial, 2) == bandwidth(D2var_banded, 2)

    # All matrix representations must contain exactly the same coefficients.
    @test D2var_full == D2var_sparse
    @test D2var_full == D2var_banded
    @test D2var_full == D2var_banded_new

    # The different matrix representations sum the same coefficients in a
    # different order, so we only get agreement up to round-off errors.
    atol = 20_000 * eps(T)
    mul!(dest1, D2var_serial, u)
    mul!(dest2, D2var_full, u)
    @test all(i -> isapprox(dest1[i], dest2[i]; atol), eachindex(u))
    mul!(dest2, D2var_sparse, u)
    @test all(i -> isapprox(dest1[i], dest2[i]; atol), eachindex(u))
    mul!(dest2, D2var_banded, u)
    @test all(i -> isapprox(dest1[i], dest2[i]; atol), eachindex(u))
end
