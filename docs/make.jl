using Documenter
import Pkg
using SummationByPartsOperators

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(SummationByPartsOperators,
                    :DocTestSetup, :(using SummationByPartsOperators); recursive = true)

# Copy some files from the top level directory to the docs and modify them
# as necessary
open(joinpath(@__DIR__, "src", "license.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/ranocha/SummationByPartsOperators.jl/blob/main/LICENSE.md"
            ```
            """)
    # Write the modified contents
    println(io, "# License")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/ranocha/SummationByPartsOperators.jl/blob/main/CODE_OF_CONDUCT.md"
            ```
            """)
    # Write the modified contents
    println(io, "# [Code of Conduct](@id code-of-conduct)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/ranocha/SummationByPartsOperators.jl/blob/main/CONTRIBUTING.md"
            ```
            """)
    # Write the modified contents
    println(io, "# Contributing")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        println(io, "> ", line)
    end
end

# Make documentation
makedocs(modules = [SummationByPartsOperators],
         sitename = "SummationByPartsOperators.jl",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                                  canonical = "https://ranocha.github.io/SummationByPartsOperators.jl/stable"),
         # Explicitly specify documentation structure
         pages = [
             "Home" => "index.md",
             "Introduction" => "introduction.md",
             "Tutorials" => [
                 "tutorials/basic_interface.md",
                 "tutorials/constant_linear_advection.md",
                 "tutorials/advection_diffusion.md",
                 "tutorials/variable_linear_advection.md",
                 "tutorials/wave_equation.md",
                 "tutorials/kdv.md",
                 "tutorials/twodimensional_linear_advection.md"
             ],
             "Automatic differentiation (AD)" => "ad.md",
             "Applications & references" => "applications.md",
             "Benchmarks" => "benchmarks.md",
             "API reference" => "api_reference.md",
             "Contributing" => "contributing.md",
             "Code of conduct" => "code_of_conduct.md",
             "License" => "license.md"
         ])

deploydocs(repo = "github.com/ranocha/SummationByPartsOperators.jl",
           devbranch = "main",
           # Only push previews if all the relevant environment variables are non-empty.
           push_preview = all(!isempty, (get(ENV, "GITHUB_TOKEN", ""), get(ENV, "DOCUMENTER_KEY", ""))))
