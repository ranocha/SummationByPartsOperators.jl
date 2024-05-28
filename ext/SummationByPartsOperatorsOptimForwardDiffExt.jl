module SummationByPartsOperatorsOptimForwardDiffExt

if isdefined(Base, :get_extension)
    using Optim: Options, LBFGS, optimize, minimizer
    using ForwardDiff: ForwardDiff
else
    using ..Optim: LBFGS, optimize, minimizer
    using ..ForwardDiff: ForwardDiff
end

using SummationByPartsOperators: SummationByPartsOperators, GlaubitzNordströmÖffner2023
using LinearAlgebra: Diagonal, diag, norm, cond

function vandermonde_matrix(functions, nodes)
    N = length(nodes)
    K = length(functions)
    V = zeros(eltype(nodes), N, K)
    for i in 1:N
        for j in 1:K
            V[i, j] = functions[j](nodes[i])
        end
    end
    return V
end

function create_S(sigma, N)
    S = zeros(eltype(sigma), N, N)
    k = 1
    for i in 1:N
        for j in (i + 1):N
            S[i, j] = sigma[k]
            k += 1
        end
    end
    return S - S'
end

function create_P(rho, x_length)
    sig(x) = 1 / (1 + exp(-x))
    P = Diagonal(sig.(rho))
    P *= x_length / sum(P)
    return P
end

function SummationByPartsOperators.construct_function_space_operator(basis_functions, x_min, x_max, nodes, ::GlaubitzNordströmÖffner2023)
    K = length(basis_functions)
    N = length(nodes)
    basis_functions_derivatives = [x -> ForwardDiff.derivative(basis_functions[i], x) for i in 1:K]
    # TODO: Orthogonalize basis w.r.t. H1 inner product via Gram-Schmidt
    V = vandermonde_matrix(basis_functions, nodes)
    V_x = vandermonde_matrix(basis_functions_derivatives, nodes)
    W = [V; -V_x]
    # display(cond(W))

    B = zeros((N, N))
    B[1, 1] = -1.0
    B[N, N] = 1.0

    R = B * V / 2
    x_length = x_max - x_min
    p = (W, R, x_length)
    function optimization_function(x, p)
        W, R, x_length = p
        (N, _) = size(R)
        L = Integer(N*(N - 1)/2)
        sigma = x[1:L]
        rho = x[(L + 1):end]
        S = create_S(sigma, N)
        P = create_P(rho, x_length)
        X = [S P]
        # TODO: Which norm to use? XW + R is a NxK matrix. This uses the Frobenius norm.
        return norm(X * W + R)^2
    end
    L = Integer(N*(N - 1)/2)
    x0 = zeros(L + N)
    result = optimize(x -> optimization_function(x, p), x0, LBFGS(), Options(g_tol = 1e-14, iterations = 10000), autodiff = :forward)
    # display(result)
    x = minimizer(result)
    sigma = x[1:L]
    rho = x[(L + 1):end]
    S = create_S(sigma, N)
    P = create_P(rho, x_length)
    weights = diag(P)
    Q = S + B/2
    D = inv(P) * Q
    return weights, D
end

end # module
