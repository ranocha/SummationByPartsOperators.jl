module SummationByPartsOperatorsOptimForwardDiffExt

using Optim: Optim, Options, LBFGS, optimize, minimizer
using ForwardDiff: ForwardDiff

using SummationByPartsOperators: SummationByPartsOperators, GlaubitzNordströmÖffner2023,
                                 MatrixDerivativeOperator
using LinearAlgebra: Diagonal, LowerTriangular, dot, diag, norm, mul!
using SparseArrays: spzeros

function SummationByPartsOperators.function_space_operator(basis_functions,
                                                           nodes::Vector{T},
                                                           source::SourceOfCoefficients;
                                                           derivative_order = 1,
                                                           accuracy_order = 0,
                                                           opt_alg = LBFGS(),
                                                           options = Options(g_tol = 1e-14,
                                                                             iterations = 10000),
                                                           verbose = false) where {T,
                                                                                   SourceOfCoefficients
                                                                                   }
    if derivative_order != 1
        throw(ArgumentError("Derivative order $derivative_order not implemented."))
    end
    sort!(nodes)
    weights, D = construct_function_space_operator(basis_functions, nodes, source;
                                                   opt_alg = opt_alg, options = options,
                                                   verbose = verbose)
    return MatrixDerivativeOperator(first(nodes), last(nodes), nodes, weights, D,
                                    accuracy_order, source)
end

function inner_H1(f, g, f_derivative, g_derivative, nodes)
    return sum(f.(nodes) .* g.(nodes) + f_derivative.(nodes) .* g_derivative.(nodes))
end
norm_H1(f, f_derivative, nodes) = sqrt(inner_H1(f, f, f_derivative, f_derivative, nodes))

function call_orthonormal_basis_function(A, basis_functions, k, x)
    sum([basis_functions[i](x) * A[k, i] for i in 1:k])
end

# This will orthonormalize the basis functions using the Gram-Schmidt process to reduce the condition
# number of the Vandermonde matrix. The matrix A transfers the old basis functions to the new orthonormalized by
# g(x) = A * f(x), where f(x) is the vector of old basis functions and g(x) is the vector of the new orthonormalized
# basis functions. Analogously, we have g'(x) = A * f'(x).
function orthonormalize_gram_schmidt(basis_functions, basis_functions_derivatives, nodes)
    K = length(basis_functions)

    A = LowerTriangular(zeros(eltype(nodes), K, K))

    basis_functions_orthonormalized = Vector{Function}(undef, K)
    basis_functions_orthonormalized_derivatives = Vector{Function}(undef, K)

    for k in 1:K
        A[k, k] = 1
        for j in 1:(k - 1)
            g(x) = call_orthonormal_basis_function(A, basis_functions, j, x)
            function g_derivative(x)
                call_orthonormal_basis_function(A, basis_functions_derivatives, j, x)
            end
            inner_product = inner_H1(basis_functions[k], g, basis_functions_derivatives[k],
                                     g_derivative, nodes)
            norm_squared = inner_H1(g, g, g_derivative, g_derivative, nodes)
            A[k, :] = A[k, :] - inner_product / norm_squared * A[j, :]
        end

        basis_functions_orthonormalized[k] = x -> call_orthonormal_basis_function(A,
                                                                                  basis_functions,
                                                                                  k, x)
        basis_functions_orthonormalized_derivatives[k] = x -> call_orthonormal_basis_function(A,
                                                                                              basis_functions_derivatives,
                                                                                              k,
                                                                                              x)
        # Normalization
        r = norm_H1(basis_functions_orthonormalized[k],
                    basis_functions_orthonormalized_derivatives[k], nodes)
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
    set_S!(S, sigma, N)
    return S
end

function set_S!(S, sigma, N)
    k = 1
    for i in 1:N
        for j in (i + 1):N
            S[i, j] = sigma[k]
            S[j, i] = -sigma[k]
            k += 1
        end
    end
end

sig(x) = 1 / (1 + exp(-x))
sig_deriv(x) = sig(x) * (1 - sig(x))

function create_P(rho, x_length)
    P = Diagonal(sig.(rho))
    P *= x_length / sum(P)
    return P
end

function construct_function_space_operator(basis_functions, nodes,
                                           ::GlaubitzNordströmÖffner2023;
                                           opt_alg = LBFGS(),
                                           options = Options(g_tol = 1e-14,
                                                             iterations = 10000),
                                           verbose = false)
    K = length(basis_functions)
    N = length(nodes)
    L = div(N * (N - 1), 2)
    basis_functions_derivatives = [x -> ForwardDiff.derivative(basis_functions[i], x)
                                   for i in 1:K]
    basis_functions_orthonormalized, basis_functions_orthonormalized_derivatives = orthonormalize_gram_schmidt(basis_functions,
                                                                                                               basis_functions_derivatives,
                                                                                                               nodes)
    V = vandermonde_matrix(basis_functions_orthonormalized, nodes)
    V_x = vandermonde_matrix(basis_functions_orthonormalized_derivatives, nodes)
    # Here, W satisfies W'*W = I
    # W = [V; -V_x]

    B = spzeros(eltype(nodes), N, N)
    B[1, 1] = -1
    B[N, N] = 1

    R = B * V / 2
    x_length = last(nodes) - first(nodes)
    S = zeros(eltype(nodes), N, N)
    A = zeros(eltype(nodes), N, K)
    SV = zeros(eltype(nodes), N, K)
    PV_x = zeros(eltype(nodes), N, K)
    daij_dsigmak = zeros(eltype(nodes), N, K, L)
    daij_drhok = zeros(eltype(nodes), N, K, N)
    p = (V, V_x, R, x_length, S, A, SV, PV_x, daij_dsigmak, daij_drhok)

    x0 = zeros(L + N)
    fg!(F, G, x) = optimization_function_and_grad!(F, G, x, p)
    result = optimize(Optim.NLSolversBase.only_fg!(fg!), x0, opt_alg, options)
    verbose && display(result)

    x = minimizer(result)
    sigma = x[1:L]
    rho = x[(L + 1):end]
    S = create_S(sigma, N)
    P = create_P(rho, x_length)
    weights = diag(P)
    Q = S + B / 2
    D = inv(P) * Q
    return weights, D
end

@views function optimization_function_and_grad!(F, G, x, p)
    V, V_x, R, x_length, S, A, SV, PV_x, daij_dsigmak, daij_drhok = p
    (N, _) = size(R)
    L = div(N * (N - 1), 2)
    sigma = x[1:L]
    rho = x[(L + 1):end]
    set_S!(S, sigma, N)
    P = create_P(rho, x_length)
    mul!(SV, S, V)
    mul!(PV_x, P, V_x)
    @. A = SV - PV_x + R
    if !isnothing(G)
        fill!(daij_dsigmak, zero(eltype(daij_dsigmak)))
        for k in axes(daij_dsigmak, 3)
            for j in axes(daij_dsigmak, 2)
                for i in axes(daij_dsigmak, 1)
                    l_tilde = k + i - N * (i - 1) + div(i * (i - 1), 2)
                    # same as above, but needs more type conversions
                    # l_tilde = Int(k + i - (i - 1) * (N - i/2))
                    if i + 1 <= l_tilde <= N
                        daij_dsigmak[i, j, k] += V[l_tilde, j]
                    else
                        C = N^2 - 3 * N + 2 * i - 2 * k + 1 / 4
                        if C >= 0
                            D = sqrt(C)
                            D_plus_one_half = D + 0.5
                            D_plus_one_half_trunc = trunc(D_plus_one_half)
                            if D_plus_one_half == D_plus_one_half_trunc
                                int_D_plus_one_half = trunc(Int, D_plus_one_half_trunc)
                                l_hat = N - int_D_plus_one_half
                                if 1 <= l_hat <= i - 1
                                    daij_dsigmak[i, j, k] -= V[l_hat, j]
                                end
                            end
                        end
                    end
                end
            end
        end
        sig_rho = sig.(rho)
        sig_deriv_rho = sig_deriv.(rho)
        sum_sig_rho = sum(sig_rho)
        for k in axes(daij_drhok, 3)
            for j in axes(daij_drhok, 2)
                for i in axes(daij_drhok, 1)
                    factor1 = x_length * V_x[i, j] / sum_sig_rho^2
                    factor = factor1 * sig_deriv_rho[k]
                    if k == i
                        daij_drhok[i, j, k] = -factor * (sum_sig_rho - sig_rho[k])
                    else
                        daij_drhok[i, j, k] = factor * sig_rho[i]
                    end
                end
            end
        end
        for k in axes(daij_dsigmak, 3)
            G[k] = 2 * dot(daij_dsigmak[:, :, k], A)
        end
        for k in axes(daij_drhok, 3)
            G[L + k] = 2 * dot(daij_drhok[:, :, k], A)
        end
    end
    if !isnothing(F)
        return norm(A)^2
    end
end

end # module
