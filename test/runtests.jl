using Base.Test

tic()

@time @testset "Periodic Operators" begin include("periodic_operators_test.jl") end
@time @testset "Non-Periodic Operators" begin include("SBP_operators_test.jl") end
@time @testset "Dissipation Operators" begin include("dissipation_operators_test.jl") end
@time @testset "Variable COefficient Operators" begin include("var_coef_operators_test.jl") end
@time @testset "Banded Matrices" begin include("banded_matrices_test.jl") end
@time @testset "Fourier Operators" begin include("fourier_operators_test.jl") end

toc()
