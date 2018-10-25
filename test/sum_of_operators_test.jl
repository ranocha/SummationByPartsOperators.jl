using Test
using LinearAlgebra
using SummationByPartsOperators

for T in (Float32, Float64)
    xmin = zero(T)
    xmax = one(T)
    N = 51
    D₁ = periodic_derivative_operator(1, 4, xmin, xmax, N)
    D₂ = periodic_derivative_operator(2, 4, xmin, xmax, N)
    D₃ = periodic_derivative_operator(3, 4, xmin, xmax, N)

    x = grid(D₁)
    u = sinpi.(x)
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))
    dest3 = fill(zero(eltype(u)), length(u))

    sum_12 = D₁ + D₂
    sum_123 = D₁ + D₂ + D₃
    sum_12_3 = sum_12 + D₃
    sum_3_12 = D₃ + sum_12

    # Compare mul! with β=0 and mul! without β.
    mul!(dest1, sum_12, u, one(T), zero(T))
    mul!(dest2, sum_12, u, one(T))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))

    # Test commutativity.
    mul!(dest1, sum_123, u, T(2))
    mul!(dest2, sum_12_3, u, T(2))
    mul!(dest3, sum_3_12, u, T(2))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    @test all(i->dest1[i] ≈ dest3[i], eachindex(u))
end
