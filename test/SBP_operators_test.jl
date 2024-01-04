using Test
using LinearAlgebra
using SummationByPartsOperators

D_test_list = (MattssonNordström2004(), MattssonSvärdNordström2004(),
                MattssonSvärdShoeybi2008(), Mattsson2014(),
                MattssonAlmquistCarpenter2014Extended(),
                MattssonAlmquistCarpenter2014Optimal(),
                MattssonAlmquistVanDerWeide2018Minimal(),
                MattssonAlmquistVanDerWeide2018Accurate(),
                DienerDorbandSchnetterTiglio2007(),
                SharanBradyLivescu2022(1//2, 2//3))

# Accuracy tests of first derivative operators.
for source in D_test_list, T in (Float32,Float64)
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    der_order = 1

    acc_order = 2
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
            summary(IOContext(devnull, :compact=>compact), D)
        end
        @test real(D) == T
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 1000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->res[i] ≈ x0[i], eachindex(res))
        # only interior
        mul!(res, D, x2)
        @test all(i->res[i] ≈ 2*x1[i], inner_indices)
        mul!(res, D, x3)
        @test any(i->!(res[i] ≈ 3*x2[i]), inner_indices)
        # integration
        k=0; @test integrate(x0, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=1; @test integrate(x1, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        @test integrate(identity, x1, D) ≈ (xmax^2 - xmin^2) / 2
        # boundary derivative
        @test derivative_left( D, x1, Val{0}()) ≈ x1[begin]
        @test derivative_right(D, x1, Val{0}()) ≈ x1[end]
    end

    acc_order = 4
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 1000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->res[i] ≈ x0[i], eachindex(res))
        if source != Mattsson2014()
            mul!(res, D, x2)
            @test all(i->res[i] ≈ 2*x1[i], eachindex(res))
        end
        # only interior
        mul!(res, D, x3)
        @test all(i->res[i] ≈ 3*x2[i], inner_indices)
        mul!(res, D, x4)
        @test all(i->res[i] ≈ 4*x3[i], inner_indices)
        mul!(res, D, x5)
        @test any(i->!(res[i] ≈ 5*x4[i]), inner_indices)
        # integration
        k=0; @test integrate(x0, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=1; @test integrate(x1, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=2; @test integrate(x2, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=3; @test integrate(x3, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        @test integrate(identity, x1, D) ≈ (xmax^2 - xmin^2) / 2
        # boundary derivative
        @test derivative_left( D, x1, Val{0}()) ≈ x1[begin]
        @test derivative_right(D, x1, Val{0}()) ≈ x1[end]
    end

    acc_order = 6
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 1000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->res[i] ≈ x0[i], eachindex(res))
        mul!(res, D, x2)
        @test all(i->res[i] ≈ 2*x1[i], eachindex(res))
        if source != Mattsson2014()
            mul!(res, D, x3)
            @test all(i->res[i] ≈ 3*x2[i], eachindex(res))
        end
        # only interior
        mul!(res, D, x4)
        @test all(i->res[i] ≈ 4*x3[i], inner_indices)
        mul!(res, D, x5)
        @test all(i->res[i] ≈ 5*x4[i], inner_indices)
        mul!(res, D, x6)
        @test all(i->res[i] ≈ 6*x5[i], inner_indices)
        # integration
        k=0; @test integrate(x0, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=1; @test integrate(x1, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=2; @test integrate(x2, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=3; @test integrate(x3, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=4; @test integrate(x4, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=5; @test integrate(x5, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=6; @test integrate(x6, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=7; @test integrate(x7, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        @test integrate(identity, x1, D) ≈ (xmax^2 - xmin^2) / 2
        # boundary derivative
        @test derivative_left( D, x1, Val{0}()) ≈ x1[begin]
        @test derivative_right(D, x1, Val{0}()) ≈ x1[end]
    end

    acc_order = 8
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 16000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->isapprox(res[i], x0[i], rtol=28000*eps(T)), eachindex(res))
        mul!(res, D, x2)
        @test all(i->isapprox(res[i], 2*x1[i], rtol=18000*eps(T)), eachindex(res))
        mul!(res, D, x3)
        @test all(i->isapprox(res[i], 3*x2[i], rtol=18000*eps(T)), eachindex(res))
        mul!(res, D, x4)
        @test all(i->res[i] ≈ 4*x3[i], inner_indices)
        # only interior
        mul!(res, D, x5)
        @test all(i->res[i] ≈ 5*x4[i], inner_indices)
        mul!(res, D, x6)
        @test all(i->isapprox(res[i], 6*x5[i], rtol=8000*eps(T)), inner_indices)
        # integration
        k=0; @test integrate(x0, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=1; @test integrate(x1, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=2; @test integrate(x2, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=3; @test integrate(x3, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=4; @test integrate(x4, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=5; @test integrate(x5, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=6; @test integrate(x6, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        k=7; @test integrate(x7, D) ≈ (xmax^(k+1)-(xmin)^(k+1))/(k+1)
        @test integrate(identity, x1, D) ≈ (xmax^2 - xmin^2) / 2
        # boundary derivative
        @test derivative_left( D, x1, Val{0}()) ≈ x1[begin]
        @test derivative_right(D, x1, Val{0}()) ≈ x1[end]
    end

    @test_throws Union{MethodError,ArgumentError} derivative_operator(source, der_order, 16, xmin, xmax, N)
end

# Accuracy tests of second derivative operators.
for source in D_test_list, T in (Float32,Float64)
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    der_order = 2

    acc_order = 2
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 5000*eps(T), eachindex(res))
        # only interior
        mul!(res, D, x2)
        @test all(i->res[i] ≈ 2*x0[i], inner_indices)
        mul!(res, D, x3)
        @test all(i->res[i] ≈ 6*x1[i], inner_indices)
        mul!(res, D, x4)
        @test any(i->!(res[i] ≈ 12*x2[i]), inner_indices)
        # boundary derivative
        @test abs(derivative_left( D, x0, Val{1}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
    end

    acc_order = 4
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 2000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 10000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->isapprox(res[i], 2*x0[i], rtol=7000*eps(T)), eachindex(res))
        # only interior
        mul!(res, D, x3)
        @test all(i->res[i] ≈ 6*x1[i], inner_indices)
        mul!(res, D, x4)
        @test all(i->res[i] ≈ 12*x2[i], inner_indices)
        mul!(res, D, x5)
        @test all(i->res[i] ≈ 20*x3[i], inner_indices)
        mul!(res, D, x6)
        @test any(i->!(res[i] ≈ 30*x4[i]), inner_indices)
        # boundary derivative
        @test abs(derivative_left( D, x0, Val{1}())) < 10*eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < 10*eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
    end

    acc_order = 6
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 10000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 20000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->isapprox(res[i], 2*x0[i], rtol=15000*eps(res[i])), eachindex(res))
        mul!(res, D, x3)
        @test all(i->isapprox(res[i], 6*x1[i], rtol=5000*eps(res[i])), eachindex(res))
        # only interior
        mul!(res, D, x4)
        @test all(i->res[i] ≈ 12*x2[i], inner_indices)
        mul!(res, D, x5)
        @test all(i->res[i] ≈ 20*x3[i], inner_indices)
        mul!(res, D, x6)
        @test all(i->res[i] ≈ 30*x4[i], inner_indices)
        mul!(res, D, x7)
        @test all(i->res[i] ≈ 42*x5[i], inner_indices)
        # boundary derivative
        @test abs(derivative_left( D, x0, Val{1}())) < 40*eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < 40*eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
        @test derivative_left( D, x4, Val{1}()) ≈ 4xmin^3
        @test derivative_right(D, x4, Val{1}()) ≈ 4xmax^3
    end

    acc_order = 8
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3
        x8 = x4 .* x4
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 10000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 22000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->isapprox(res[i], 2*x0[i], rtol=40000*eps(T)), eachindex(res))
        mul!(res, D, x3)
        @test all(i->isapprox(res[i], 6*x1[i], rtol=20000*eps(T)), eachindex(res))
        mul!(res, D, x4)
        @test all(i->res[i] ≈ 12*x2[i], inner_indices)
        # only interior
        mul!(res, D, x5)
        @test all(i->res[i] ≈ 20*x3[i], inner_indices)
        mul!(res, D, x6)
        @test all(i->isapprox(res[i], 30*x4[i], rtol=10000*eps(T)), inner_indices)
        mul!(res, D, x7)
        @test all(i->res[i] ≈ 42*x5[i], inner_indices)
        mul!(res, D, x8)
        @test all(i->isapprox(res[i], 56*x6[i], rtol=20000*eps(T)), inner_indices)
        # boundary derivative
        @test abs(derivative_left( D, x0, Val{1}())) < 40*eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < 40*eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
        @test derivative_left( D, x4, Val{1}()) ≈ 4xmin^3
        @test derivative_right(D, x4, Val{1}()) ≈ 4xmax^3
        @test derivative_left( D, x5, Val{1}()) ≈ 5xmin^4
        @test derivative_right(D, x5, Val{1}()) ≈ 5xmax^4
    end

    @test_throws Union{MethodError,ArgumentError} derivative_operator(source, der_order, 16, xmin, xmax, N)
end

# Accuracy tests of third derivative operators.
@testset "third-derivative operators" for source in D_test_list, T in (Float32,Float64)
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    der_order = 3

    acc_order = 2
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 50_000*eps(T), eachindex(res))
        # only interior
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 160_000*eps(T), inner_indices)
        mul!(res, D, x3)
        @test all(i->abs(res[i] - 6*x0[i]) < 340_000*eps(T), inner_indices)
        # boundary: first derivative
        @test abs(derivative_left( D, x0, Val{1}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        # boundary: second derivative
        @test abs(derivative_left( D, x0, Val{2}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{2}())) < eps(T)
        @test abs(derivative_left( D, x1, Val{2}())) < 600 * eps(T)
        @test abs(derivative_right(D, x1, Val{2}())) < 600 * eps(T)
        @test derivative_left( D, x2, Val{2}()) ≈ 2
        @test derivative_right(D, x2, Val{2}()) ≈ 2
    end

    acc_order = 4
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x3 .* x2
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 100_000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 500_000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 500_000*eps(T), eachindex(res))
        # only interior
        mul!(res, D, x3)
        @test all(i->isapprox(res[i], 6*x0[i], rtol=100_000*eps(T)), inner_indices)
        mul!(res, D, x4)
        @test all(i->isapprox(res[i], 24*x1[i], rtol=50_000*eps(T)), inner_indices)
        mul!(res, D, x5)
        @test all(i->isapprox(res[i], 60*x2[i], rtol=50_000*eps(T)), inner_indices)
        # boundary: first derivative
        @test abs(derivative_left( D, x0, Val{1}())) < 10*eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < 10*eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ 1
        @test derivative_right(D, x1, Val{1}()) ≈ 1
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
        # boundary: second derivative
        @test abs(derivative_left( D, x0, Val{2}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{2}())) < eps(T)
        @test abs(derivative_left( D, x1, Val{2}())) < 5_000*eps(T)
        @test abs(derivative_right(D, x1, Val{2}())) < 5_000*eps(T)
        @test derivative_left( D, x2, Val{2}()) ≈ 2
        @test derivative_right(D, x2, Val{2}()) ≈ 2
        @test derivative_left( D, x3, Val{2}()) ≈ 6xmin
        @test derivative_right(D, x3, Val{2}()) ≈ 6xmax
    end

    acc_order = 6
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D) |> collect
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x3 .* x2
        x6 = x3 .* x3
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 1_000_000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 5_000_000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 5_000_000*eps(T), eachindex(res))
        mul!(res, D, x3)
        @test all(i->isapprox(res[i], 6*x0[i], rtol=5_000_000*eps(T)), eachindex(res))
        # only interior
        mul!(res, D, x3)
        @test all(i->isapprox(res[i], 6*x0[i], rtol=150_000*eps(T)), inner_indices)
        mul!(res, D, x4)
        @test all(i->isapprox(res[i], 24*x1[i], rtol=50_000*eps(T)), inner_indices)
        mul!(res, D, x5)
        @test all(i->isapprox(res[i], 60*x2[i], rtol=50_000*eps(T)), inner_indices)
        mul!(res, D, x6)
        @test all(i->isapprox(res[i], 120*x3[i], rtol=50_000*eps(T)), inner_indices)
        # boundary: first derivative
        @test abs(derivative_left( D, x0, Val{1}())) < 50*eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < 50*eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ 1
        @test derivative_right(D, x1, Val{1}()) ≈ 1
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
        # boundary: second derivative
        @test abs(derivative_left( D, x0, Val{2}())) < 5_000*eps(T)
        @test abs(derivative_right(D, x0, Val{2}())) < 5_000*eps(T)
        @test abs(derivative_left( D, x1, Val{2}())) < 10_000*eps(T)
        @test abs(derivative_right(D, x1, Val{2}())) < 11_000*eps(T)
        @test isapprox(derivative_left( D, x2, Val{2}()), 2, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x2, Val{2}()), 2, atol=20_000*eps(T))
        @test isapprox(derivative_left( D, x3, Val{2}()), 6xmin, atol=50_000*eps(T))
        @test isapprox(derivative_right(D, x3, Val{2}()), 6xmax, atol=80_000*eps(T))
    end

    @test_throws Union{MethodError,ArgumentError} derivative_operator(source, der_order, 16, xmin, xmax, N)
end

# Accuracy tests of fourth derivative operators.
@testset "fourth-derivative operators" for source in D_test_list, T in (Float32,Float64)
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    der_order = 4

    acc_order = 2
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 10 * eps(T) / D.Δx^4, eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 100 * eps(T) / D.Δx^4, eachindex(res))
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 100 * eps(T) / D.Δx^4, eachindex(res))
        # only interior
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 100 * eps(T) / D.Δx^4, inner_indices)
        mul!(res, D, x3)
        @test all(i->abs(res[i]) < 100 * eps(T) / D.Δx^4, inner_indices)
        # boundary: first derivative
        @test abs(derivative_left( D, x0, Val{1}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{1}())) < eps(T)
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        # boundary: second derivative
        @test abs(derivative_left( D, x0, Val{2}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{2}())) < eps(T)
        @test abs(derivative_left( D, x1, Val{2}())) < 600 * eps(T)
        @test abs(derivative_right(D, x1, Val{2}())) < 600 * eps(T)
        @test derivative_left( D, x2, Val{2}()) ≈ 2 * one(T)
        @test derivative_right(D, x2, Val{2}()) ≈ 2 * one(T)
        # boundary: third derivative
        @test abs(derivative_left( D, x0, Val{3}())) < 10 * eps(T) / D.Δx^3
        @test abs(derivative_right(D, x0, Val{3}())) < 10 * eps(T) / D.Δx^3
        @test abs(derivative_left( D, x1, Val{3}())) < 10 * eps(T) / D.Δx^3
        @test abs(derivative_right(D, x1, Val{3}())) < 10 * eps(T) / D.Δx^3
        @test abs(derivative_left( D, x2, Val{3}())) < 10 * eps(T) / D.Δx^3
        @test abs(derivative_right(D, x2, Val{3}())) < 10 * eps(T) / D.Δx^3
    end

    acc_order = 4
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 50_000_000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 50_000_000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 50_000_000*eps(T), eachindex(res))
        # only interior
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 50_000_000*eps(T), inner_indices)
        mul!(res, D, x3)
        @test all(i->abs(res[i]) < 120_000_000*eps(T), inner_indices)
        mul!(res, D, x4)
        @test all(i->isapprox(res[i], 24*x0[i], rtol=50_000_000*eps(T)), inner_indices)
        # boundary: first derivative
        @test isapprox(derivative_left( D, x0, Val{1}()), 0, atol=10*eps(T))
        @test isapprox(derivative_right(D, x0, Val{1}()), 0, atol=10*eps(T))
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
        @test isapprox(derivative_left( D, x4, Val{1}()), 4xmin^3, atol=1.e-3)
        @test isapprox(derivative_right(D, x4, Val{1}()), 4xmax^3, atol=1.e-3)
        # boundary: second derivative
        @test isapprox(derivative_left( D, x0, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x0, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_left( D, x1, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x1, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_left( D, x2, Val{2}()), 2, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x2, Val{2}()), 2, atol=10_000*eps(T))
        @test isapprox(derivative_left( D, x3, Val{2}()), 6xmin, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x3, Val{2}()), 6xmax, atol=50_000*eps(T))
        # boundary: third derivative
        @test abs(derivative_left( D, x0, Val{3}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{3}())) < eps(T)
        @test abs(derivative_left( D, x1, Val{3}())) < 100_000*eps(T)
        @test abs(derivative_right(D, x1, Val{3}())) < 200_000*eps(T)
        @test abs(derivative_left( D, x2, Val{3}())) < 100_000*eps(T)
        @test abs(derivative_right(D, x2, Val{3}())) < 100_000*eps(T)
        @test isapprox(derivative_left( D, x3, Val{3}()), 6, atol=150_000*eps(T))
        @test isapprox(derivative_right(D, x3, Val{3}()), 6, atol=1_000_000*eps(T))
    end

    acc_order = 6
    D = try
        derivative_operator(source, der_order, acc_order, xmin, xmax, N)
    catch
        nothing
    end
    if D !== nothing
        D = derivative_operator(source, der_order, acc_order, xmin, xmax, N)
        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D)
            show(IOContext(devnull, :compact=>compact), D.coefficients)
        end
        x1 = grid(D)
        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        res = fill(zero(eltype(x0)), length(x0))
        inner_indices = length(D.coefficients.left_boundary)+1:length(res)-length(D.coefficients.left_boundary)-1
        @test derivative_order(D) == der_order
        @test accuracy_order(D)   == acc_order
        @test issymmetric(D) == false
        # interior and boundary
        mul!(res, D, x0)
        @test all(i->abs(res[i]) < 50_000_000*eps(T), eachindex(res))
        mul!(res, D, x1)
        @test all(i->abs(res[i]) < 50_000_000*eps(T), eachindex(res))
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 100_000_000*eps(T), eachindex(res))
        # only interior
        mul!(res, D, x2)
        @test all(i->abs(res[i]) < 100_000_000*eps(T), inner_indices)
        mul!(res, D, x3)
        @test all(i->abs(res[i]) < 500_000_000*eps(T), inner_indices)
        mul!(res, D, x4)
        @test all(i->isapprox(res[i], 24*x0[i], rtol=50_000_000*eps(T)), inner_indices)
        # boundary: first derivative
        @test isapprox(derivative_left( D, x0, Val{1}()), 0, atol=1000*eps(T))
        @test isapprox(derivative_right(D, x0, Val{1}()), 0, atol=1000*eps(T))
        @test derivative_left( D, x1, Val{1}()) ≈ one(T)
        @test derivative_right(D, x1, Val{1}()) ≈ one(T)
        @test derivative_left( D, x2, Val{1}()) ≈ 2xmin
        @test derivative_right(D, x2, Val{1}()) ≈ 2xmax
        @test derivative_left( D, x3, Val{1}()) ≈ 3xmin^2
        @test derivative_right(D, x3, Val{1}()) ≈ 3xmax^2
        @test isapprox(derivative_left( D, x4, Val{1}()), 4xmin^3, atol=1.e-3)
        @test isapprox(derivative_right(D, x4, Val{1}()), 4xmax^3, atol=1.e-3)
        # boundary: second derivative
        @test isapprox(derivative_left( D, x0, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x0, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_left( D, x1, Val{2}()), 0, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x1, Val{2}()), 0, atol=11_000*eps(T))
        @test isapprox(derivative_left( D, x2, Val{2}()), 2, atol=10_000*eps(T))
        @test isapprox(derivative_right(D, x2, Val{2}()), 2, atol=20_000*eps(T))
        @test isapprox(derivative_left( D, x3, Val{2}()), 6xmin, atol=20_000*eps(T))
        @test isapprox(derivative_right(D, x3, Val{2}()), 6xmax, atol=75_000*eps(T))
        # boundary: third derivative
        @test abs(derivative_left( D, x0, Val{3}())) < eps(T)
        @test abs(derivative_right(D, x0, Val{3}())) < eps(T)
        @test abs(derivative_left( D, x1, Val{3}())) < 500_000*eps(T)
        @test abs(derivative_right(D, x1, Val{3}())) < 500_000*eps(T)
        @test abs(derivative_left( D, x2, Val{3}())) < 500_000*eps(T)
        @test abs(derivative_right(D, x2, Val{3}())) < 1_000_000*eps(T)
        @test isapprox(derivative_left( D, x3, Val{3}()), 6, atol=900_000*eps(T))
        @test isapprox(derivative_right(D, x3, Val{3}()), 6, atol=5_000_000*eps(T))
    end

    @test_throws Union{MethodError,ArgumentError} derivative_operator(source, der_order, 16, xmin, xmax, N)
end


# Compare mul! with β=0 and mul! without β.
for T in (Float32, Float64), acc_order in (2,4,6,8)
    xmin = zero(T)
    xmax = 5*one(T)
    N = 51
    source = MattssonSvärdShoeybi2008()
    D_serial = derivative_operator(source, 1, acc_order, xmin, xmax, N)
    D_threads= derivative_operator(source, 1, acc_order, xmin, xmax, N, ThreadedMode())
    D_safe   = derivative_operator(source, 1, acc_order, xmin, xmax, N, SafeMode())
    D_full   = Matrix(D_serial)
    x = grid(D_serial)
    u = x.^5
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    mul!(dest1, D_serial, u, one(T), zero(T))
    mul!(dest2, D_serial, u, one(T))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest1, D_safe, u, one(T), zero(T))
    mul!(dest2, D_safe, u, one(T))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest1, D_threads, u, one(T), zero(T))
    mul!(dest2, D_threads, u, one(T))
    @test all(i->dest1[i] ≈ dest2[i], eachindex(u))
    dest3 = D_serial*u
    @test all(i->dest1[i] ≈ dest3[i], eachindex(u))
    dest3 = D_full*u
    @test all(i->dest1[i] ≈ dest3[i], eachindex(u))
    dest3 = D_safe*u
    @test all(i->dest1[i] ≈ dest3[i], eachindex(u))
end


# https://github.com/ranocha/SummationByPartsOperators.jl/issues/208
@testset "issue #208" begin
    for acc_order in (2, 4, 6), der_order in 1:4
        @test_nowarn derivative_operator(Mattsson2014(),
                                         derivative_order = der_order,
                                         accuracy_order = acc_order,
                                         xmin = -1.0, xmax = 1.0,
                                         N = 20)
    end
end
