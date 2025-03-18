using Aqua
using SummationByPartsOperators

Aqua.test_all(SummationByPartsOperators;
  ambiguities = false, # a lot of false positives from dependencies
  unbound_args = false, # TODO: a strange problem I do not understand right now
  stale_deps = (; ignore = [:PrecompileTools, :InteractiveUtils]),
)
