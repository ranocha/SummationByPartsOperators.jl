# Benchmark `mul!` with first-derivative SBP operators for a variety of element
# types. The kernels of SummationByPartsOperators.jl treat composite element
# types such as `StaticVector`s, `ForwardDiff.Dual`s, and `Complex` numbers
# differently from plain scalars, see
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
#
# The script is self-contained: it activates a temporary environment and
# installs fixed versions of all relevant packages together with the requested
# version of SummationByPartsOperators.jl. Hence, no setup is required and the
# results can be reproduced on other machines.
#
#
# ## Measuring
#
# Pass the branch (or tag or commit) of SummationByPartsOperators.jl to measure
# and the file the results shall be written to:
#
#     julia --threads=1 benchmark/benchmark_element_types.jl main results_main.csv
#     julia --threads=1 benchmark/benchmark_element_types.jl hr/dual results_hr_dual.csv
#
# Each run takes roughly 15 to 20 minutes with the default settings. The results
# contain the name of the machine, the versions of all relevant packages, and
# `versioninfo()`, so they can be copied between computers.
#
#
# ## Comparing
#
# Pass two result files to get a summary of the speedups and a merged file
# `comparison.csv` containing every single ratio:
#
#     julia benchmark/benchmark_element_types.jl --compare results_main.csv results_hr_dual.csv
#
# This does not install anything. Only compare results obtained on the same
# machine. Comparing two runs of the *same* version is a good way to get a
# feeling for the noise level of a machine; treat speedups within that range as
# noise.
#
#
# ## Options (environment variables)
#
# - `SBP_BENCHMARK_SECONDS`: time budget per measurement in seconds (default 0.2)
# - `SBP_BENCHMARK_REPEATS`: how often the whole grid is repeated; the minimum
#   over all repetitions is reported (default 2)
# - `SBP_BENCHMARK_MODES`: comma-separated subset of `FastMode`, `SafeMode`,
#   `ThreadedMode` (default `FastMode,SafeMode`)
# - `SBP_BENCHMARK_QUICK`: set to `1` for a much smaller grid to test the setup

const REPOSITORY = "https://github.com/ranocha/SummationByPartsOperators.jl"

# Fixed versions of the packages that influence the measurements so that
# different versions of SummationByPartsOperators.jl are compared on equal
# footing - also across machines and at later points in time.
# `LoopVectorization` is not used by this script directly but determines how the
# kernels of `FastMode` are vectorized.
const DEPENDENCIES = ["BenchmarkTools" => v"1.8.0",
    "ForwardDiff" => v"1.4.6",
    "StaticArrays" => v"1.9.22",
    "LoopVectorization" => v"0.12.174"]

const USAGE = """
              Usage:
                julia --threads=1 $(PROGRAM_FILE) BRANCH RESULTS.csv
                julia $(PROGRAM_FILE) --compare OLD.csv NEW.csv

              where BRANCH is a branch, tag, or commit of
              $(REPOSITORY), e.g., `main` or `hr/dual`."""

function parse_arguments(args)
    if length(args) == 3 && args[1] == "--compare"
        return (:compare, args[2], args[3])
    elseif length(args) == 2 && !startswith(args[1], "-")
        return (:run, args[1], args[2])
    else
        return (:usage, "", "")
    end
end

# ------------------------------------------------------------------------------
# Comparing two result files. This part must not use any package beyond the
# standard library so that comparing results does not require an installation.
# ------------------------------------------------------------------------------

using Printf: @printf
using Statistics: median

function read_results(path)
    results = Dict{NTuple{6, String}, Float64}()
    order = NTuple{6, String}[]
    for line in eachline(path)
        (isempty(line) || startswith(line, '#') || startswith(line, "elementtype")) &&
            continue
        fields = split(line, ',')
        length(fields) == 9 || continue
        key = ntuple(i -> String(fields[i]), 6)
        haskey(results, key) || push!(order, key)
        # keep the minimum over all repetitions
        results[key] = min(get(results, key, Inf), parse(Float64, fields[7]))
    end
    return results, order
end

function print_metadata(name, path)
    println("# ", name, ": ", path)
    for line in eachline(path)
        startswith(line, '#') || break
        any(startswith(line, prefix)
            for prefix in ("# host", "# SummationByPartsOperators.jl branch",
                           "# package version", "# julia threads",
                           "# benchmark seconds")) || continue
        println("  ", line)
    end
    return nothing
end

function compare_results(path_old, path_new)
    old, order = read_results(path_old)
    new, _ = read_results(path_new)
    common = [key for key in order if haskey(new, key)]
    isempty(common) && error("no common configurations in $path_old and $path_new")

    print_metadata("old", path_old)
    print_metadata("new", path_new)
    println("\nSpeedup = time(old) / time(new); > 1 means the new version is faster.")
    println("Aggregated over the numbers of nodes, accuracy orders, (non)periodic ",
            "operators,\nand 3-arg/5-arg `mul!`.\n")

    labels = unique(key[1] for key in common)
    modes = unique(key[2] for key in common)
    @printf("%-18s", "element type")
    foreach(mode -> @printf("%31s", mode), modes)
    println()
    @printf("%-18s", "")
    for _ in modes
        @printf("%8s%9s%8s%6s", "min", "median", "max", "<0.9")
    end
    println()
    for label in labels
        @printf("%-18s", label)
        for mode in modes
            ratios = [old[key] / new[key]
                      for key in common if key[1] == label && key[2] == mode]
            isempty(ratios) && continue
            @printf("%8.2f%9.2f%8.2f%4d/%-2d", minimum(ratios), median(ratios),
                    maximum(ratios), count(<(0.9), ratios), length(ratios))
        end
        println()
    end

    open("comparison.csv", "w") do io
        println(io,
                "elementtype,mode,nnodes,accuracy_order,operator,call," *
                "time_old_ns,time_new_ns,speedup")
        for key in common
            println(io, join((key..., old[key], new[key], old[key] / new[key]), ','))
        end
    end
    println("\nWrote comparison.csv with all ", length(common), " ratios.")
    return nothing
end

# ------------------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------------------

const COMMAND, FIRST_ARGUMENT, SECOND_ARGUMENT = parse_arguments(ARGS)

if COMMAND === :usage
    println(USAGE)
    exit(1)
elseif COMMAND === :compare
    compare_results(FIRST_ARGUMENT, SECOND_ARGUMENT)
else
    using Pkg
    Pkg.activate(; temp = true)
    packages = [Pkg.PackageSpec(name = name, version = version)
                for (name, version) in DEPENDENCIES]
    push!(packages, Pkg.PackageSpec(url = REPOSITORY, rev = FIRST_ARGUMENT))
    Pkg.add(packages)

    # The measurements are in a separate file since they use macros of
    # BenchmarkTools.jl, which is only available after the installation above.
    include("measure_element_types.jl")
    run_benchmarks(FIRST_ARGUMENT, SECOND_ARGUMENT)
end
