
"""
    GlaubitzNordströmÖffner2023()

Function space SBP (FSBP) operators given in
- Glaubitz, Nordström, Öffner (2023)
  Summation-by-parts operators for general function spaces.
  SIAM Journal on Numerical Analysis 61, 2, pp. 733-754.

See also
- Glaubitz, Nordström, Öffner (2024)
  An optimization-based construction procedure for function space based
  summation-by-parts operators on arbitrary grids.
  arXiv, arXiv:2405.08770v1.

See [`function_space_operator`](@ref).
"""
struct GlaubitzNordströmÖffner2023 <: SourceOfCoefficients end

function Base.show(io::IO, source::GlaubitzNordströmÖffner2023)
  if get(io, :compact, false)
    summary(io, source)
  else
      print(io,
          "Glaubitz, Nordström, Öffner (2023) \n",
          "  Summation-by-parts operators for general function spaces \n",
          "  SIAM Journal on Numerical Analysis 61, 2, pp. 733-754. \n",
          "See also \n",
          "  Glaubitz, Nordström, Öffner (2024) \n",
          "  An optimization-based construction procedure for function \n",
          "    space based summation-by-parts operators on arbitrary grids \n",
          "  arXiv, arXiv:2405.08770v1.")
  end
end

"""
    function_space_operator(basis_functions, x_min, x_max, nodes, source; accuracy_order = 0)

Construct an operator that represents a derivative operator in a function space spanned by
the `basis_functions`, which is an iterable of functions. The operator is constructed on the
interval `[x_min, x_max]` with the nodes `nodes`. The `accuracy_order` is the order of the
accuracy of the operator, which can optionally be passed, but does not have any effect on the
operator.
The operator that is returned follows the general interface. Currently, it is wrapped in a
[`MatrixDerivativeOperator`](@ref), but this might change in the future.

See also [`GlaubitzNordströmÖffner2023`](@ref).
"""
function function_space_operator(basis_functions,
                                 x_min::T, x_max::T, nodes::Vector{T},
                                 source::SourceOfCoefficients;
                                 accuracy_order = 0) where {T, SourceOfCoefficients}

  weights, D = construct_function_space_operator(basis_functions, x_min, x_max, nodes, source)
  return MatrixDerivativeOperator(x_min, x_max, nodes, weights, D, accuracy_order, source)
end

# This function is extended in the package extension SummationByPartsOperatorsOptimExt
function construct_function_space_operator end
