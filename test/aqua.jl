using Aqua
using SummationByPartsOperators

Aqua.test_all(SummationByPartsOperators;
  ambiguities = false, # a lot of false positives from dependencies
  unbound_args = false, # TODO: a strange problem I do not understand right now
  stale_deps = (; ignore = [:PrecompileTools]),
  # We would like to test the Project.toml formatting but there are some
  # CI issues, see https://github.com/ranocha/BSeries.jl/pull/119
  project_toml_formatting = false,
)

# Project.toml formatting only on newer versions of Julia
if VERSION >= v"1.9"
    Aqua.test_project_toml_formatting(BSeries)
end
