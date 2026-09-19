# The measurements of `benchmark_element_types.jl`. This file is included by
# that script after it has installed the packages used below; it is not meant to
# be run on its own.

using BenchmarkTools
using ForwardDiff: ForwardDiff
using InteractiveUtils: versioninfo
using LinearAlgebra: mul!
using Random: Xoshiro
using StaticArrays: SVector
using SummationByPartsOperators

const QUICK = get(ENV, "SBP_BENCHMARK_QUICK", "0") == "1"
const NUMBERS_OF_NODES = QUICK ? [30] : [30, 100, 1000]
const ACCURACY_ORDERS = QUICK ? [4] : [2, 4, 6]
const REPEATS = parse(Int, get(ENV, "SBP_BENCHMARK_REPEATS", QUICK ? "1" : "2"))
const SECONDS = parse(Float64, get(ENV, "SBP_BENCHMARK_SECONDS", "0.2"))
const CSV_HEADER = "elementtype,mode,nnodes,accuracy_order,operator,call," *
                   "time_min_ns,time_median_ns,allocations"

# The element types to benchmark together with the real number type their values
# are made of. The labels must not contain commas since they are written to a
# CSV file.
function element_type_cases()
    tag = :SummationByPartsOperatorsBenchmark
    Dual = ForwardDiff.Dual
    return [("Float32", Float32, Float32),
            ("Float64", Float64, Float64),
            ("ComplexF32", ComplexF32, Float32),
            ("ComplexF64", ComplexF64, Float64),
            ("SVector{1}", SVector{1, Float64}, Float64),
            ("SVector{2}", SVector{2, Float64}, Float64),
            ("SVector{3}", SVector{3, Float64}, Float64),
            ("SVector{4}", SVector{4, Float64}, Float64),
            ("SVector{5}", SVector{5, Float64}, Float64),
            ("SVector{8}", SVector{8, Float64}, Float64),
            ("Dual{1}", Dual{tag, Float64, 1}, Float64),
            ("Dual{2}", Dual{tag, Float64, 2}, Float64),
            ("Dual{3}", Dual{tag, Float64, 3}, Float64),
            ("Dual{4}", Dual{tag, Float64, 4}, Float64),
            ("Dual{5}", Dual{tag, Float64, 5}, Float64),
            ("Dual{8}", Dual{tag, Float64, 8}, Float64),
            ("Dual{12}", Dual{tag, Float64, 12}, Float64),
            ("Dual32{4}", Dual{tag, Float32, 4}, Float32),
            ("Dual{Dual{2}}{2}", Dual{:outer, Dual{tag, Float64, 2}, 2}, Float64)]
end

function execution_modes()
    available = Dict("FastMode" => FastMode(), "SafeMode" => SafeMode(),
                     "ThreadedMode" => ThreadedMode())
    names = split(get(ENV, "SBP_BENCHMARK_MODES", "FastMode,SafeMode"), ',')
    return [available[strip(name)] for name in names]
end

# All element types considered here are isbits values consisting of
# `num_components` numbers of type `R` stored contiguously. Thus, random values
# can be created - and all components can be compared - without knowing anything
# else about them.
num_components(::Type{T}, ::Type{R}) where {T, R} = div(sizeof(T), sizeof(R))

function random_vector(::Type{T}, ::Type{R}, N, rng) where {T, R}
    values = rand(rng, R, num_components(T, R) * N)
    return collect(reinterpret(T, values))
end

function components(u::AbstractVector{T}, ::Type{R}) where {T, R}
    return num_components(T, R) == 1 ? reinterpret(R, u) : vec(reinterpret(reshape, R, u))
end

function make_operators(::Type{R}, mode, N, accuracy_order) where {R}
    D = derivative_operator(MattssonNordström2004(); derivative_order = 1,
                            accuracy_order = accuracy_order, xmin = zero(R),
                            xmax = one(R), N = N, mode = mode)
    D_periodic = periodic_derivative_operator(derivative_order = 1,
                                              accuracy_order = accuracy_order,
                                              xmin = zero(R), xmax = one(R), N = N,
                                              mode = mode)
    return ("nonperiodic" => D, "periodic" => D_periodic)
end

"""
    check_correctness(::Type{T}, ::Type{R}, mode, N, accuracy_order)

Compare the results of `mul!` to the dense matrix representation of the
operators. This guards against benchmarking a broken kernel. All components of
the values are compared since, e.g., `isapprox` of vectors of `Dual`s ignores
the partial derivatives.
"""
function check_correctness(::Type{T}, ::Type{R}, mode, N, accuracy_order) where {T, R}
    rng = Xoshiro(1234)
    u = random_vector(T, R, N, rng)
    du = similar(u)
    α = convert(R, 2)
    β = convert(R, 3)
    rtol = sqrt(eps(R))

    for (_, D) in make_operators(R, mode, N, accuracy_order)
        A = Matrix(D)
        mul!(du, D, u)
        isapprox(components(du, R), components(A * u, R); rtol = rtol) || return false

        dest = copy(u)
        mul!(dest, D, u, α, β)
        isapprox(components(dest, R), components(α * (A * u) + β * u, R); rtol = rtol) ||
            return false
    end
    return true
end

function benchmark_case(io, label, ::Type{T}, ::Type{R}, mode, N,
                        accuracy_order) where {T, R}
    rng = Xoshiro(1234)
    u = random_vector(T, R, N, rng)
    du = similar(u)
    α = convert(R, 2)
    β = convert(R, 3)
    mode_name = string(nameof(typeof(mode)))

    for (operator_name, D) in make_operators(R, mode, N, accuracy_order)
        fill!(du, zero(T))
        mul!(du, D, u) # compile before measuring
        trial3 = @benchmark mul!($du, $D, $u)

        fill!(du, zero(T))
        mul!(du, D, u, α, β)
        fill!(du, zero(T))
        trial5 = @benchmark mul!($du, $D, $u, $α, $β)

        for (call, trial) in (("3arg", trial3), ("5arg", trial5))
            best = minimum(trial)
            println(io,
                    join((label, mode_name, N, accuracy_order, operator_name, call,
                          BenchmarkTools.time(best), BenchmarkTools.time(median(trial)),
                          BenchmarkTools.allocs(best)), ','))
        end
    end
    flush(io)
    return nothing
end

function installed_versions()
    versions = Pair{String, String}[]
    for (_, package) in Pkg.dependencies()
        package.version === nothing && continue
        (package.is_direct_dep || package.name == "LoopVectorization") || continue
        push!(versions, package.name => string(package.version))
    end
    return sort!(versions; by = first)
end

function write_header(io, rev)
    println(io, "# SummationByPartsOperators.jl element type benchmark")
    println(io, "# date: ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    println(io, "# host: ", Base.gethostname())
    println(io, "# SummationByPartsOperators.jl branch: ", rev)
    for (name, version) in installed_versions()
        println(io, "# package version: ", name, " ", version)
    end
    println(io, "# julia threads: ", Threads.nthreads())
    println(io, "# benchmark seconds: ", SECONDS, ", repeats: ", REPEATS)
    buffer = IOBuffer()
    versioninfo(buffer)
    for line in eachline(IOBuffer(take!(buffer)))
        println(io, "# ", line)
    end
    println(io, CSV_HEADER)
    flush(io)
    return nothing
end

function run_benchmarks(rev, path)
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = SECONDS
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 10_000
    cases = element_type_cases()
    modes = execution_modes()

    print("Checking correctness ")
    for (label, T, R) in cases, mode in modes
        check_correctness(T, R, mode, 30, 4) ||
            error("$label is not computed correctly using $mode")
        print(".")
    end
    println(" ok")

    total = REPEATS * length(cases) * length(modes) * length(NUMBERS_OF_NODES) *
            length(ACCURACY_ORDERS)
    done = 0
    start = time()
    open(path, "w") do io
        write_header(io, rev)
        for _ in 1:REPEATS
            for (label, T, R) in cases, mode in modes, N in NUMBERS_OF_NODES,
                accuracy_order in ACCURACY_ORDERS

                benchmark_case(io, label, T, R, mode, N, accuracy_order)
                done += 1
                @printf("\r%4d/%4d configurations, %6.1f s elapsed   ", done, total,
                        time()-start)
            end
        end
    end
    println("\nWrote ", path)
    return nothing
end
