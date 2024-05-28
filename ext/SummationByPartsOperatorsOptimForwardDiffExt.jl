module SummationByPartsOperatorsOptimForwardDiffExt

if isdefined(Base, :get_extension)
    using Optim: Options, LBFGS, optimize, minimizer
    using ForwardDiff: ForwardDiff
else
    using ..Optim: LBFGS, optimize, minimizer
    using ..ForwardDiff: Options, ForwardDiff
end

using SummationByPartsOperators: SummationByPartsOperators, GlaubitzNordströmÖffner2023
using LinearAlgebra: Diagonal, LowerTriangular, diag, norm, cond

function inner_H1(f, g, f_derivative, g_derivative, nodes)
    return sum(f.(nodes) .* g.(nodes) + f_derivative.(nodes) .* g_derivative.(nodes))
end
norm_H1(f, f_derivative, nodes) = sqrt(inner_H1(f, f, f_derivative, f_derivative, nodes))

call_orthonormal_basis_function(A, basis_functions, k, x) = sum([basis_functions[i](x)*A[k, i] for i in 1:k])

# This will orthonormalize the basis functions using the Gram-Schmidt process to reduce the condition
# number of the Vandermonde matrix. The matrix A transfers the old basis functions to the new orthonormalized by
# g(x) = A * f(x), where f(x) is the vector of old basis functions and g(x) is the vector of the new orthonormalized
# basis functions. Analogously, we have g'(x) = A * f'(x).
function orthonormalize_gram_schmidt(basis_functions, basis_functions_derivatives, nodes)
    K = length(basis_functions)

    A = LowerTriangular(zeros(eltype(nodes), K, K))

    basis_functions_orthonormalized = Vector{Function}(undef, K)
    basis_functions_orthonormalized_derivatives = Vector{Function}(undef, K)

    for k = 1:K
        A[k, k] = 1.0
        for j = 1:k-1
            g = x -> call_orthonormal_basis_function(A, basis_functions, j, x)
            g_derivative = x -> call_orthonormal_basis_function(A, basis_functions_derivatives, j, x)
            inner_product = inner_H1(basis_functions[k], g, basis_functions_derivatives[k], g_derivative, nodes)
            norm_squared = inner_H1(g, g, g_derivative, g_derivative, nodes)
            A[k, :] = A[k, :] - inner_product/norm_squared * A[j, :]
        end

        basis_functions_orthonormalized[k] = x -> call_orthonormal_basis_function(A, basis_functions, k, x)
        basis_functions_orthonormalized_derivatives[k] = x -> call_orthonormal_basis_function(A, basis_functions_derivatives, k, x)
        # Normalization
        r = norm_H1(basis_functions_orthonormalized[k], basis_functions_orthonormalized_derivatives[k], nodes)
        A[k, :] = A[k, :] / r
    end
    return basis_functions_orthonormalized, basis_functions_orthonormalized_derivatives
end

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

sig(x) = 1 / (1 + exp(-x))

function create_P(rho, x_length)
    P = Diagonal(sig.(rho))
    P *= x_length / sum(P)
    return P
end

function SummationByPartsOperators.construct_function_space_operator(basis_functions, x_min, x_max, nodes,
                                                                     ::GlaubitzNordströmÖffner2023)
    K = length(basis_functions)
    N = length(nodes)
    basis_functions_derivatives = [x -> ForwardDiff.derivative(basis_functions[i], x) for i in 1:K]
    basis_functions_orthonormalized, basis_functions_orthonormalized_derivatives = orthonormalize_gram_schmidt(basis_functions,
                                                                                                               basis_functions_derivatives,
                                                                                                               nodes)
    V = vandermonde_matrix(basis_functions_orthonormalized, nodes)
    V_x = vandermonde_matrix(basis_functions_orthonormalized_derivatives, nodes)
    # Here, W satisfies W'*W = I
    W = [V; -V_x]

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
    result = optimize(x -> optimization_function(x, p), x0, LBFGS(), Options(g_tol = 1e-14, iterations = 10000),
                      autodiff = :forward)
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
