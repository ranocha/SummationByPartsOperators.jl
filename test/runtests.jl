using Test

const SBP_TEST = get(ENV, "SBP_TEST", "all")

@time @testset "SummationByPartsOperators.jl tests" begin
    @time if SBP_TEST == "all" || SBP_TEST == "part1"
        @time @testset "Periodic Operators" begin include("periodic_operators_test.jl") end
        @time @testset "Non-Periodic Operators" begin include("SBP_operators_test.jl") end
        @time @testset "Dissipation Operators" begin include("dissipation_operators_test.jl") end
        @time @testset "Variable Coefficient Operators" begin include("var_coef_operators_test.jl") end
        @time @testset "Banded Matrices" begin include("banded_matrices_test.jl") end
        @time @testset "Fourier Operators" begin include("fourier_operators_test.jl") end
        @time @testset "Fourier Operators 2D" begin include("fourier_operators_2d_test.jl") end
        @time @testset "Legendre Operators" begin include("legendre_operators_test.jl") end
        @time @testset "Linear Combinations of Operators" begin include("linear_combinations_of_operators_test.jl") end
        @time @testset "Upwind Operators" begin include("upwind_operators_test.jl") end
        @time @testset "Special Matrix Types" begin include("special_matrix_types.jl") end
        @time @testset "Aqua" begin include("aqua.jl") end
    end

    @time if SBP_TEST == "all" || SBP_TEST == "part2"
        @time @testset "Coupling" begin include("coupling_test.jl") end
        @time @testset "Conservation Laws" begin
            include("conservation_laws/burgers_test.jl")
            include("conservation_laws/cubic_test.jl")
            include("conservation_laws/variable_linear_advection_test.jl")
            include("conservation_laws/quartic_nonconvex_test.jl")
        end
        @time @testset "Second Order Equations" begin
            include("second_order_eqs/wave_eq_test.jl")
        end
    end
end
