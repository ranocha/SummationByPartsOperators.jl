using Test

@time @testset "Periodic Operators" begin include("periodic_operators_test.jl") end
@time @testset "Non-Periodic Operators" begin include("SBP_operators_test.jl") end
@time @testset "Dissipation Operators" begin include("dissipation_operators_test.jl") end
@time @testset "Variable Coefficient Operators" begin include("var_coef_operators_test.jl") end
@time @testset "Banded Matrices" begin include("banded_matrices_test.jl") end
@time @testset "Fourier Operators" begin include("fourier_operators_test.jl") end
@time @testset "Legendre Operators" begin include("legendre_operators_test.jl") end
@time @testset "Sum of Operators" begin include("sum_of_operators_test.jl") end
@time @testset "Conservation Laws" begin
    include("conservation_laws/burgers_test.jl")
    include("conservation_laws/cubic_test.jl")
    include("conservation_laws/variable_linear_advection_test.jl")
end
@time @testset "Second Order Equations" begin
    include("second_order_eqs/wave_eq_test.jl")
end
