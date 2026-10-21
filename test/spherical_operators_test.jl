using LinearAlgebra
using SparseArrays
using Test

@testset "Spherical SBP operators" begin
    source = VretinarisSchnetter2026()

    for grid_kind in (:collocated, :staggered), order in (4, 6)
        operators = spherical_operators(source;
                                        accuracy_order = order,
                                        N = 20,
                                        R = 20 // 1,
                                        p = 2,
                                        grid = grid_kind,
                                        mode = SafeMode())

        r = grid(operators)
        expected_grid = grid_kind === :collocated ?
                        Rational{BigInt}.(0:20) :
                        Rational{BigInt}.((1:2:39) .* 20 // 39)
        @test r == expected_grid
        @test SummationByPartsOperators.accuracy_order(operators) == order
        @test source_of_coefficients(operators) === source
        G = gradient_operator(operators)
        Godd = odd_derivative_operator(operators)
        D = divergence_operator(operators)
        @test size(G) == (length(r), length(r))
        @test size(D) == (length(r), length(r))
        @test eltype(G) == eltype(r)
        @test grid(G) === r
        @test accuracy_order(D) == order
        @test source_of_coefficients(G) === source
        @test_throws ArgumentError SphericalOperatorView(operators, operators.gradient,
                                                         :invalid)
        @test Matrix(G) == operators.gradient
        @test sparse(D) == operators.divergence
        gradient = G * ones(eltype(r), length(r))
        gradient_inplace = similar(gradient)
        mul!(gradient_inplace, G, ones(eltype(r), length(r)))
        @test gradient_inplace == gradient
        @test scalar_mass(operators) * D + G' * vector_mass(operators) ==
              boundary_matrix(operators)
        @test boundary_matrix(operators)[end, end] == 20^2
        @test gradient == zeros(eltype(r), length(r))

        divergence = D * r
        if grid_kind === :collocated
            exact_rows = order == 4 ? (1:(length(r) - 4)) : eachindex(r)
            @test divergence[exact_rows] == 3 .* ones(eltype(r), length(exact_rows))
            @test scalar_mass(operators) != vector_mass(operators)
        else
            @test divergence[1:16] == 3 .* ones(eltype(r), 16)
            @test scalar_mass(operators) == vector_mass(operators)
        end
    end
end

@testset "Spherical SBP operators with p=1" begin
    source = VretinarisSchnetter2026()

    for grid_kind in (:collocated, :staggered), order in (4, 6)
        operators = spherical_operators(source;
                                        accuracy_order = order,
                                        N = 20,
                                        R = 20 // 1,
                                        p = 1,
                                        grid = grid_kind,
                                        mode = SafeMode())
        r = grid(operators)
        G = gradient_operator(operators)
        D = divergence_operator(operators)

        @test scalar_mass(operators) * D + G' * vector_mass(operators) ==
              boundary_matrix(operators)
        @test boundary_matrix(operators)[end, end] == 20
        @test G * ones(eltype(r), length(r)) == zeros(eltype(r), length(r))
        if grid_kind === :collocated && order == 4
            Dr = D * r
            Dr3 = D * (r .^ 3)
            @test Dr[1:(end - 4)] == 2 .* ones(eltype(r), length(r) - 4)
            @test Dr3[1:5] == 4 .* r[1:5] .^ 2
        end
    end
end

@testset "Fourth-order cylindrical spherical convergence" begin
    gradient_errors = Float64[]
    divergence_errors = Float64[]
    for N in (20, 40, 80)
        operators = spherical_operators(VretinarisSchnetter2026();
                                        accuracy_order = 4,
                                        N,
                                        R = 1.0,
                                        p = 1,
                                        mode = SafeMode())
        r = grid(operators)
        G = gradient_operator(operators)
        D = divergence_operator(operators)
        interior = 10:(length(r) - 9)

        push!(gradient_errors, maximum(abs.(G * cos.(r) .+ sin.(r))[interior]))
        flux = r .* exp.(-r .^ 2)
        exact_divergence = (2 .- 2 .* r .^ 2) .* exp.(-r .^ 2)
        push!(divergence_errors, maximum(abs.(D * flux .- exact_divergence)[interior]))
    end
    @test all(gradient_errors[1:(end - 1)] ./ gradient_errors[2:end] .> 8)
    @test all(divergence_errors[1:(end - 1)] ./ divergence_errors[2:end] .> 8)
end

@testset "Spherical SBP input validation" begin
    source = VretinarisSchnetter2026()
    @test_throws ArgumentError spherical_operators(source;
                                                   accuracy_order = 2,
                                                   N = 20,
                                                   R = 1,
                                                   grid = :staggered)
    @test_throws ArgumentError spherical_operators(source;
                                                   accuracy_order = 4,
                                                   N = 0,
                                                   R = 1,
                                                   grid = :staggered)
    @test_throws ArgumentError spherical_operators(source;
                                                   accuracy_order = 4,
                                                   N = 20,
                                                   R = 1,
                                                   grid = :invalid)
    error = try
        spherical_operators(source;
                            accuracy_order = 4,
                            N = 20,
                            R = 1,
                            p = 0)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    @test occursin("usual Cartesian operators", sprint(showerror, error))
    @test occursin("derivative_operator", sprint(showerror, error))
    @test_throws ArgumentError spherical_operators(source;
                                                   accuracy_order = 4,
                                                   N = 20,
                                                   R = 1,
                                                   p = 3)
end
