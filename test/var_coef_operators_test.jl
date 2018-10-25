using Test
using LinearAlgebra, SparseArrays
using SummationByPartsOperators

test_list = (Mattsson2012(),)

# Test consistency with constant coefficient operators
for source in test_list, acc_order in (2,4,6), T in (Float32,Float64)
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101

    D2 = derivative_operator(source, 2, acc_order, xmin, xmax, N)
    D2var = try
        var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one)
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    D2var == nothing && continue

    println(devnull, D2var)
    println(devnull, D2var.coefficients)
    @test maximum(abs, Matrix(D2) - Matrix(D2var)) < 10000*eps(T)
end


# Compare mul! with β=0 and mul! without β.
for T in (Float32, Float64), acc_order in (2,4,6)
    xmin = zero(T)
    xmax = 5*one(T)
    N = 51
    source = Mattsson2012()

    D2var_serial = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one, Val{:serial}())
    D2var_threads= var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one, Val{:threads}())
    D2var_full   = Matrix(D2var_serial)
    D2var_sparse = sparse(D2var_serial)

    x = grid(D2var_serial)
    u = x.^5
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    mul!(dest1, D2var_serial, u, one(T), zero(T))
    mul!(dest2, D2var_serial, u, one(T))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest1, D2var_threads, u, one(T), zero(T))
    mul!(dest2, D2var_threads, u, one(T))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest2, D2var_serial, u)
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest2, D2var_full, u)
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
end
