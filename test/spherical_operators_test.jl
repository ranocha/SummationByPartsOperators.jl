using LinearAlgebra
using SparseArrays
using Test

@testset "Spherical SBP operators" begin
    source = VretinarisSchnetter2026()

    for grid_kind in (:collocated, :staggered), order in (2, 4, 6)
        operators = spherical_operators(source;
                                        accuracy_order = order,
                                        N = 20,
                                        R = 20 // 1,
                                        p = 2,
                                        grid = grid_kind)

        r = grid(operators)
        expected_grid = grid_kind === :collocated ?
                        Rational{BigInt}.(0:20) :
                        Rational{BigInt}.((1:2:39) .* 20 // 39)
        @test r == expected_grid
        @test SummationByPartsOperators.accuracy_order(operators) == order
        @test source_of_coefficients(operators) === source
        G = gradient_operator(operators)
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
            if order == 2
                @test scalar_mass(operators) == vector_mass(operators)
            else
                @test scalar_mass(operators) != vector_mass(operators)
            end
        else
            @test scalar_mass(operators) == vector_mass(operators)
            exact_rows = order == 2 ? eachindex(r) : 1:16
            @test divergence[exact_rows] == 3 .* ones(eltype(r), length(exact_rows))
        end
    end
end

@testset "Spherical SBP operators with p=1" begin
    source = VretinarisSchnetter2026()

    for grid_kind in (:collocated, :staggered), order in (2, 4, 6)
        operators = spherical_operators(source;
                                        accuracy_order = order,
                                        N = 20,
                                        R = 20 // 1,
                                        p = 1,
                                        grid = grid_kind)
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

@testset "Second-order spherical closure" begin
    for grid_kind in (:collocated, :staggered), p in (1, 2)
        operators = spherical_operators(VretinarisSchnetter2026();
                                        accuracy_order = 2,
                                        N = 20,
                                        R = 20 // 1,
                                        p,
                                        grid = grid_kind)
        r = grid(operators)
        S = scalar_mass(operators)
        V = vector_mass(operators)
        B = boundary_matrix(operators)
        G = gradient_operator(operators)
        Godd = operators.odd_derivative
        D = divergence_operator(operators)

        @test S == V
        @test isdiag(S)
        @test isdiag(V)
        @test all(>(zero(eltype(r))), diag(S))
        @test S * D + G' * V == B
        @test D * r == (p + 1) .* ones(eltype(r), length(r))
        @test sum(diag(S)) == last(r)^(p + 1) / (p + 1)
        if grid_kind === :collocated
            @test D[1, :] == (p + 1) .* Godd[1, :]
        end
    end
end

@testset "Second-order spherical convergence" begin
    for grid_kind in (:collocated, :staggered)
        gradient_errors = Float64[]
        divergence_errors = Float64[]
        for N in (20, 40, 80)
            operators = spherical_operators(VretinarisSchnetter2026();
                                            accuracy_order = 2,
                                            N,
                                            R = 1.0,
                                            p = 2,
                                            grid = grid_kind)
            r = grid(operators)
            G = gradient_operator(operators)
            D = divergence_operator(operators)
            interior = 2:(length(r) - 1)

            push!(gradient_errors, maximum(abs.(G * cos.(r) .+ sin.(r))[interior]))
            flux = r .* exp.(-r .^ 2)
            exact_divergence = (3 .- 2 .* r .^ 2) .* exp.(-r .^ 2)
            push!(divergence_errors,
                  maximum(abs.(D * flux .- exact_divergence)[interior]))
        end
        @test all(gradient_errors[1:(end - 1)] ./ gradient_errors[2:end] .> 3)
        @test all(divergence_errors[1:(end - 1)] ./ divergence_errors[2:end] .> 3)
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
                                        p = 1)
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

@testset "Supported Cartesian sources for spherical operators" begin
    supported = ((MattssonNordström2004(), (2, 4, 6)),
                 (MattssonSvärdNordström2004(), (2, 4, 6)),
                 (MattssonSvärdShoeybi2008(), (2, 4, 6)),
                 (Mattsson2012(), (2, 4, 6)),
                 (Mattsson2014(), (2,)),
                 (MattssonAlmquistCarpenter2014Extended(), (2,)),
                 (MattssonAlmquistCarpenter2014Optimal(), (2,)),
                 (Mattsson2017(:minus), (4, 6)),
                 (DienerDorbandSchnetterTiglio2007(), (2, 4)))

    for (cartesian_source, orders) in supported
        @test SummationByPartsOperators._spherical_supported_orders(cartesian_source) ==
              orders
        for order in orders
            operators = spherical_operators(VretinarisSchnetter2026();
                                            accuracy_order = order,
                                            N = 20,
                                            R = 20 // 1,
                                            p = 2,
                                            grid = :collocated,
                                            cartesian_source)
            @test operators.cartesian_source === cartesian_source
            @test scalar_mass(operators) * divergence_operator(operators) +
                  gradient_operator(operators)' * vector_mass(operators) ==
                  boundary_matrix(operators)
        end
    end

    unsupported = (Mattsson2017(:central),
                   Mattsson2017(:plus),
                   MattssonAlmquistVanDerWeide2018Minimal(),
                   MattssonAlmquistVanDerWeide2018Accurate(),
                   MattssonNiemeläWinters2026(:central),
                   MattssonNiemeläWinters2026(:plus),
                   MattssonNiemeläWinters2026(:minus),
                   SharanBradyLivescu2022(1 // 2, 1 // 2),
                   WilliamsDuru2024(:central),
                   WilliamsDuru2024(:plus),
                   WilliamsDuru2024(:minus))
    for cartesian_source in unsupported
        @test isempty(SummationByPartsOperators._spherical_supported_orders(cartesian_source))
        @test_throws ArgumentError spherical_operators(VretinarisSchnetter2026();
                                                       accuracy_order = 4,
                                                       N = 20,
                                                       R = 1.0,
                                                       cartesian_source)
    end
    @test_throws ArgumentError spherical_operators(VretinarisSchnetter2026();
                                                   accuracy_order = 2,
                                                   N = 20,
                                                   R = 1.0,
                                                   cartesian_source = Mattsson2017(:minus))
end

@testset "Spherical SBP output element type" begin
    source = VretinarisSchnetter2026()
    radii_and_types = ((Float32(1), Float32),
                       (1.0, Float64),
                       (big"1.0", BigFloat),
                       (1, Float64),
                       (BigInt(1), BigFloat),
                       (1 // 1, Rational{BigInt}))

    for (R, T) in radii_and_types
        operators = spherical_operators(source;
                                        accuracy_order = 2,
                                        N = 8,
                                        R)
        @test eltype(grid(operators)) == T
        @test eltype(scalar_mass(operators)) == T
        @test eltype(vector_mass(operators)) == T
        @test eltype(gradient_operator(operators)) == T
        @test eltype(divergence_operator(operators)) == T
    end
end

@testset "Spherical SBP input validation" begin
    source = VretinarisSchnetter2026()
    @test_throws ArgumentError spherical_operators(source;
                                                   accuracy_order = 3,
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
