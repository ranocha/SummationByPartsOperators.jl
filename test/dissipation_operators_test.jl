using Test
using LinearAlgebra
using SummationByPartsOperators
using OffsetArrays: OffsetArray

D_test_list = (MattssonNordström2004(), MattssonSvärdNordström2004(),
               MattssonSvärdShoeybi2008(), Mattsson2014(),
               MattssonAlmquistCarpenter2014Extended(),
               MattssonAlmquistCarpenter2014Optimal(),
               MattssonAlmquistVanDerWeide2018Minimal(),
               MattssonAlmquistVanDerWeide2018Accurate())
Di_test_list = (MattssonSvärdNordström2004(),)

# Nonperiodic: Test symmetry and eigenvalues. The dissipation operators are
# adapted to the norm and the grid of `D`.
for source_D in D_test_list, source_Di in Di_test_list, acc_order in 2:2:12,
    T in (Float32, Float64)

    xmin = -one(T)
    xmax = 2 * one(T)
    N = 101
    der_order = 1

    D = try
        derivative_operator(source_D, der_order, acc_order, xmin, xmax, N)
    catch err
        !isa(err, ArgumentError) && throw(err)
        nothing
    end
    D === nothing && continue

    @inferred mass_matrix(D)
    H = mass_matrix(D)

    for order in 2:2:8
        Di = try
            dissipation_operator(source_Di, D, order = order)
        catch err
            !isa(err, ArgumentError) && throw(err)
            nothing
        end
        Di === nothing && continue

        println(devnull, Di)
        println(devnull, Di.coefficients)
        HDi = H * Matrix(Di)
        # On grids with nonuniform spacing near the boundaries, the entries of
        # `HDi` in the boundary rows are larger, and so are the rounding errors.
        scale = grid(D) isa AbstractRange ? one(T) : norm(HDi)
        @test norm(HDi - HDi') < 10 * eps(T) * scale
        @test minimum(real, eigvals(HDi)) < 10 * eps(T)
    end
end

# Nonperiodic: Compare mul! with β=0 and mul! without β.
for T in (Float32, Float64), order in (2, 4, 6, 8)
    xmin = zero(T)
    xmax = 5 * one(T)
    N = 51
    D = derivative_operator(MattssonSvärdNordström2004(), 1, order, xmin, xmax, N)
    x = grid(D)
    u = x .^ 5
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    Di_serial = dissipation_operator(D, order = order, mode = FastMode())
    Di_threads = dissipation_operator(D, order = order, mode = ThreadedMode())
    Di_safe = dissipation_operator(D, order = order, mode = SafeMode())
    Di_full = Matrix(Di_serial)

    mul!(dest1, Di_serial, u, one(T), zero(T))
    mul!(dest2, Di_serial, u, one(T))
    @test isapprox(dest1, dest2, rtol = 120 * sqrt(eps(T)))
    mul!(dest1, Di_safe, u, one(T), zero(T))
    mul!(dest2, Di_safe, u, one(T))
    @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest1, Di_threads, u, one(T), zero(T))
    mul!(dest2, Di_threads, u, one(T))
    @test isapprox(dest1, dest2, rtol = 120 * sqrt(eps(T)))
    dest3 = Di_serial * u
    @test isapprox(dest1, dest3, rtol = 120 * sqrt(eps(T)))
    dest3 = Di_safe * u
    @test all(i -> isapprox(dest1[i], dest3[i], atol = 800000 * eps(T)), eachindex(u))
    dest3 = Di_full * u
    @test all(i -> isapprox(dest1[i], dest3[i], atol = 800000 * eps(T)), eachindex(u))
end

# Nonperiodic: Dissipation operators on grids with nonuniform spacing near the
# boundaries.
function dissipation_mul_allocations(dest, Di, u)
    mul!(dest, Di, u)
    @allocated mul!(dest, Di, u)
end
function dissipation_mul_allocations(dest, Di, u, α, β)
    mul!(dest, Di, u, α, β)
    @allocated mul!(dest, Di, u, α, β)
end

for (make_source, acc_orders) in ((T -> MattssonAlmquistCarpenter2014Optimal(),
                                   (2, 4, 6, 8)),
                                  (T -> MattssonAlmquistVanDerWeide2018Minimal(),
                                   (4, 6, 8)),
                                  (T -> MattssonAlmquistVanDerWeide2018Accurate(),
                                   (4, 6, 8)),
                                  (T -> SharanBradyLivescu2022(T(1 // 2), T(2 // 3)),
                                   (2, 4, 6))),
    acc_order in acc_orders, T in (Float32, Float64)

    source_D = make_source(T)
    xmin = -one(T)
    xmax = 2 * one(T)
    p = acc_order ÷ 2

    maxerrs = T[]
    for N in (41, 81)
        D = derivative_operator(source_D, 1, acc_order, xmin, xmax, N)
        x = grid(D)
        H = mass_matrix(D)
        Di_serial = dissipation_operator(D, mode = FastMode())
        Di_threads = dissipation_operator(D, mode = ThreadedMode())
        Di_safe = dissipation_operator(D, mode = SafeMode())
        @test grid(Di_serial) == x

        Di_full = Matrix(Di_serial)
        @test Di_full ≈ Matrix(Di_threads)
        @test Di_full ≈ Matrix(Di_safe)
        HDi = H * Di_full
        @test norm(HDi - HDi') < 10 * eps(T) * norm(HDi)
        @test maximum(eigvals(Symmetric(HDi))) < 10 * eps(T) * opnorm(HDi)

        # The dissipation operator annihilates polynomials of degree less
        # than `p`.
        tol = 10 * eps(T) * opnorm(Di_full, Inf)
        for k in 0:(p - 1)
            @test maximum(abs, Di_serial * x .^ k) < tol * maximum(abs, x)^k
        end

        u = sinpi.(x)
        dest1 = Di_serial * u
        dest2 = fill(one(T), length(u))
        mul!(dest2, Di_serial, u, 2 * one(T), 3 * one(T))
        @test dest2 ≈ 2 .* dest1 .+ 3
        @test maximum(abs, Di_safe * u - dest1) < tol
        @test maximum(abs, Di_threads * u - dest1) < tol
        @test maximum(abs, Di_full * u - dest1) < tol
        @test iszero(dissipation_mul_allocations(dest1, Di_serial, u))
        @test iszero(dissipation_mul_allocations(dest2, Di_serial, u, one(T), zero(T)))
        u_offset = OffsetArray(u, -1)
        @test_throws ArgumentError mul!(similar(u_offset), Di_serial, u_offset)
        @test_throws ArgumentError mul!(similar(u_offset), Di_serial, u_offset, 2, 3)

        push!(maxerrs, maximum(abs, dest1))
    end
    # The dissipation is of order `p` near the boundaries (and of order
    # `2p` in the interior).
    if T == Float64
        @test log2(maxerrs[1] / maxerrs[2]) > p - 0.1
    end
end

# Dense reference implementation of `-H⁻¹ D_pᵀ B D_p` with `H = Δx * W`, where
# the row of `D_p` acting on the nodes `j:(j + p)` is `p! Δx^p` times the
# `p`-th divided difference on these nodes and `B = I` for these rows.
function dissipation_matrix_reference(D, order)
    x = grid(D)
    N = length(x)
    Δx = step(x)
    p = order ÷ 2
    w = diag(mass_matrix(D)) ./ Δx
    A = zeros(eltype(x), N, N)
    for j in 1:(N - p)
        idx = j:(j + p)
        d = [factorial(p) * Δx^p / prod(x[k] - x[m] for m in idx if m != k)
             for k in idx]
        A[idx, idx] .+= d * d'
    end
    return -(A ./ w)
end

# Nonperiodic: Dissipation operators on small grids with nonuniform spacing
# near the boundaries, including grids where all rows are boundary rows.
for (source_D, acc_orders) in ((MattssonAlmquistCarpenter2014Optimal(), (2, 4, 6, 8)),
                               (MattssonAlmquistVanDerWeide2018Minimal(),
                                (4, 6, 8, 10, 12)),
                               (MattssonAlmquistVanDerWeide2018Accurate(),
                                (4, 6, 8, 10, 12)),
                               (SharanBradyLivescu2022(0.5, 2 / 3), (2, 4, 6))),
    acc_order in acc_orders, order in 2:2:8

    N_min = findfirst(N -> try
                          derivative_operator(source_D, 1, acc_order, -1.0, 2.0, N)
                          true
                      catch err
                          err isa Union{ArgumentError, DimensionMismatch} || throw(err)
                          false
                      end, 1:100)
    for N in N_min:(N_min + 2 * order)
        D = derivative_operator(source_D, 1, acc_order, -1.0, 2.0, N)
        for mode in (SafeMode(), FastMode(), ThreadedMode())
            Di = dissipation_operator(D; order, mode)
            @test grid(Di) == grid(D)
            @test Matrix(Di) ≈ dissipation_matrix_reference(D, order)
        end
    end
end

# Periodic: Test symmetry and eigenvalues
for acc_order in 2:2:8, T in (Float32, Float64)
    xmin = -one(T)
    xmax = 2 * one(T)
    N = 100
    der_order = 1

    D = periodic_derivative_operator(der_order, acc_order, xmin, xmax, N)
    for order in 2:2:8
        Di = dissipation_operator(D, order = order)
        println(devnull, Di)
        HDi = Matrix(Di)
        @test norm(HDi - HDi') < 10 * eps(T)
        @test minimum(real, eigvals(HDi)) < 10 * eps(T)
    end
end

# Periodic: Compare mul! with β=0 and mul! without β.
for T in (Float32, Float64), order in (2, 4, 6, 8)
    xmin = zero(T)
    xmax = 5 * one(T)
    N = 50
    D = periodic_derivative_operator(1, order, xmin, xmax, N)
    x = grid(D)
    u = x .^ 5
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    Di_serial = dissipation_operator(D, order = order, mode = FastMode())
    Di_threads = dissipation_operator(D, order = order, mode = ThreadedMode())
    Di_safe = dissipation_operator(D, order = order, mode = SafeMode())
    Di_full = Matrix(Di_serial)

    mul!(dest1, Di_serial, u, one(T), zero(T))
    mul!(dest2, Di_serial, u, one(T))
    @test maximum(abs, dest1 - dest2) < u[end] * N * eps(T)
    mul!(dest1, Di_threads, u, one(T), zero(T))
    mul!(dest2, Di_threads, u, one(T))
    @test maximum(abs, dest1 - dest2) < u[end] * N * eps(T)
    dest3 = Di_serial * u
    @test maximum(abs, dest1 - dest3) < u[end] * N * eps(T)
    dest3 = Di_full * u
    @test maximum(abs, dest1 - dest3) < 5 * u[end] * N * eps(T)
    dest3 = Di_safe * u
    @test maximum(abs, dest1 - dest3) < 5 * u[end] * N * eps(T)
end

# Compare periodic and nonperiodic.
for T in (Float32, Float64), order in (2, 4, 6, 8)
    xmin = zero(T)
    xmax = 5 * one(T)
    N = 50
    D = derivative_operator(MattssonSvärdNordström2004(), 1, order, xmin, xmax, N + 1)
    Dp = periodic_derivative_operator(1, order, xmin, xmax, N)

    Di = dissipation_operator(D)
    Dip = dissipation_operator(Dp)

    @test norm(Matrix(Di)[15:35, 15:35] - Matrix(Dip)[15:35, 15:35]) < eps(T)
end
