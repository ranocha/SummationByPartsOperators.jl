module SummationByPartsOperatorsOptimExt

if isdefined(Base, :get_extension)
    using Optim: LBFGS, optimize
else
    using ..Optim: LBFGS, optimize
end

using SummationByPartsOperators: SummationByPartsOperators, GlaubitzNordströmÖffner2023

function SummationByPartsOperators.construct_function_space_operator(basis_functions, x_min, x_max, nodes, ::GlaubitzNordströmÖffner2023)
    # TODO: implement the construction of the function space operator by solving an optimization problem
    N = length(nodes)
    weights = ones(N)
    D = zeros(N, N)
    return weights, D
end

end # module
