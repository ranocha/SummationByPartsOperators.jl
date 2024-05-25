
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
    FunctionSpaceOperator(basis_functions, x_min, x_max, nodes, source)

Construct a `MatrixDerivativeOperator` that represents a derivative operator in a function space.
"""
function FunctionSpaceOperator(basis_functions,
                               x_min::T, x_max::T, nodes::Vector{T},
                               source::SourceOfCoefficients) where {T, SourceOfCoefficients}

  weights, D = construct_function_space_operator(basis_functions, x_min, x_max, nodes, source)
  accuracy_order = length(basis_functions) # TODO: this might needs to be adjusted
  return MatrixDerivativeOperator(x_min, x_max, nodes, weights, D, accuracy_order, source)
end

# This function is extended in the package extension SummationByPartsOperatorsOptimExt
function construct_function_space_operator end
