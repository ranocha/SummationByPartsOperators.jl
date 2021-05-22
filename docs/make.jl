using Documenter
import Pkg
using SummationByPartsOperators

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(SummationByPartsOperators,
  :DocTestSetup, :(using SummationByPartsOperators); recursive=true)

open(joinpath(@__DIR__, "src", "license.md"), "w") do io
  println(io, "# License")
  println(io, "")
  for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
    line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
    println(io, "> ", line)
  end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
  println(io, "# Contributing")
  println(io, "")
  for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
    line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
    println(io, "> ", line)
  end
end

# Make documentation
makedocs(
  modules = [SummationByPartsOperators],
  sitename="SummationByPartsOperators.jl",
  format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical = "https://ranocha.github.io/SummationByPartsOperators.jl/stable"
  ),
  # Explicitly specify documentation structure
  pages = [
    "Home" => "index.md",
    "Introduction" => "introduction.md",
    "Tutorials" => [
      "tutorials/constant_linear_advection.md",
      "tutorials/advection_diffusion.md",
      "tutorials/variable_linear_advection.md",
      # "tutorials/wave_equation.md", TODO
    ],
    "Applications" => "applications.md",
    "Benchmarks" => "benchmarks.md",
    "API reference" => "api_reference.md",
    "Contributing" => "contributing.md",
    "License" => "license.md"
  ],
  strict = true # to make the GitHub action fail when doctests fail
)

deploydocs(
  repo = "github.com/ranocha/SummationByPartsOperators.jl",
  devbranch = "main",
  push_preview = true
)
