using Base.Test

tic()

@time @testset "Periodic Operators" begin include("periodic_operators_test.jl") end

toc()
