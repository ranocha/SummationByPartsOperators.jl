
"""
    FourierDerivativeOperator{T<:Real, GridCompute, GridEvaluate, RFFT, IRFFT}

A derivative operator on a periodic grid with scalar type `T` computing the
first derivative using a spectral Fourier expansion via real discrete Fourier
transforms.
"""
struct FourierDerivativeOperator{T<:Real, GridCompute, GridEvaluate, RFFT, IRFFT} <: AbstractDerivativeOperator{T}
    jac::T
    grid_compute::GridCompute   # N-1 nodes, including the left and excluding the right boundary
    grid_evaluate::GridEvaluate #  N  nodes, including both boundaries
    tmp::Vector{Complex{T}}
    rfft_plan::RFFT
    irfft_plan::IRFFT

    function FourierDerivativeOperator(jac::T, grid_compute::GridCompute, grid_evaluate::GridEvaluate,
                                        tmp::Vector{Complex{T}}, rfft_plan::RFFT, irfft_plan::IRFFT) where {T<:Real, GridCompute, GridEvaluate, RFFT, IRFFT}
        @argcheck length(irfft_plan) == length(tmp)
        @argcheck length(irfft_plan) == (length(rfft_plan)÷2)+1
        @argcheck length(grid_compute) == length(rfft_plan)
        @argcheck length(grid_compute) == length(grid_evaluate)-1
        @argcheck first(grid_compute) == first(grid_evaluate)
        @argcheck step(grid_compute) ≈ step(grid_evaluate)
        @argcheck last(grid_compute) < last(grid_evaluate)

        new{T, GridCompute, GridEvaluate, RFFT, IRFFT}(jac, grid_compute, grid_evaluate, tmp, rfft_plan, irfft_plan)
    end
end

"""
    FourierDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}

Construct the `FourierDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` Fourier modes.
"""
function FourierDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}
    @argcheck N >= 1

    jac = 2*T(π) / (xmax - xmin)
    grid_evaluate = linspace(xmin, xmax, N+1)
    grid_compute = linspace(xmin, grid_evaluate[end-1], N)
    u = zero.(grid_compute)
    rfft_plan = plan_rfft(u)
    uhat = rfft_plan*u
    irfft_plan = plan_irfft(uhat, N)

    FourierDerivativeOperator(jac, grid_compute, grid_evaluate, uhat, rfft_plan, irfft_plan)
end

function fourier_derivative_operator(xmin::T, xmax::T, N::Int) where {T<:Real}
    FourierDerivativeOperator(xmin, xmax, N)
end

derivative_order(D::FourierDerivativeOperator) = 1
Base.issymmetric(D::FourierDerivativeOperator) = false
grid(D::FourierDerivativeOperator) = D.grid_compute

function Base.show(io::IO, D::FourierDerivativeOperator{T}) where {T}
    grid = D.grid_evaluate
    print(io, "Periodic 1st derivative Fourier operator {T=", T, "} \n")
    print(io, "on a grid in [", first(grid), ", ", last(grid),
                "] using ", length(grid)-1, " modes. \n")
end


"""
    compute_coefficients(u, D::FourierDerivativeOperator)

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D`.
"""
function compute_coefficients(u, D::FourierDerivativeOperator)
    grid = D.grid_compute
    xmin = first(grid)
    xmax = last(grid)
    uval = Array{typeof(u((xmin+xmax)/2))}(length(grid))
    compute_coefficients!(uval, u, D)
end

"""
    compute_coefficients!(uval::AbstractVector, u, D::FourierDerivativeOperator

Compute the nodal values of the function `u` at the grid associated to the
derivative operator `D` and stores the result in `uval`.
"""
function compute_coefficients!(uval::AbstractVector, u, D::FourierDerivativeOperator)
    uval .= u.(D.grid_compute)
end


"""
    evaluate_coefficients(u, D::FourierDerivativeOperator)

Evaluates the nodal coefficients `u` at a uniform grid associated to the Fourier
derivative operator `D` including both endpoints. Returns `xplot, uplot`,
where `xplot` contains the equally spaced nodes and `uplot` the corresponding
values of `u`.
"""
function evaluate_coefficients(u, D::FourierDerivativeOperator)
    grid = D.grid_evaluate
    xplot = Array{eltype(grid)}(length(grid))
    uplot = Array{eltype(u)}(length(grid))

    evaluate_coefficients!(xplot, uplot, u, D)
end

"""
    evaluate_coefficients!(xplot, uplot, u, D::FourierDerivativeOperator)

Evaluates the nodal coefficients `u` at a uniform grid associated to the Fourier
derivative operator `D` including both endpoints and store the result in `xplot,
uplot`. Returns `xplot, uplot`, where `xplot` contains the equally spaced nodes
and `uplot` the corresponding values of `u`.
"""
function evaluate_coefficients!(xplot, uplot, u, D::FourierDerivativeOperator)
    @argcheck length(uplot) == length(xplot)
    @argcheck length(uplot) == length(D.grid_evaluate)

    xplot .= D.grid_evaluate
    uplot[1:end-1] = u
    uplot[end] = uplot[1]

    xplot, uplot
end


function Base.A_mul_B!(dest::AbstractVector{T}, D::FourierDerivativeOperator,
                        u::AbstractVector{T}) where {T}
    @unpack jac, tmp, rfft_plan, irfft_plan = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    A_mul_B!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        tmp[j] *= (j-1)*im * jac
    end
    @inbounds tmp[end] = 0
    A_mul_B!(dest, irfft_plan, tmp)

    nothing
end


"""
    fourier_derivative_matrix(N, xmin::Real=0.0, xmax::Real=2π)

Compute the Fourier derivative matrix with respect to the corresponding nodal
basis, see Kopriva (2009) Implementing Spectral Methods for PDEs, Algorithm 18.
"""
function fourier_derivative_matrix(N, xmin::Real=0.0, xmax::Real=2π)
    T = promote_type(typeof(xmin), typeof(xmax))
    jac_2 = T(π) / (xmax - xmin)
    D = Array{T}(N, N)
    @inbounds for j in 1:N, i in 1:N
        j == i && continue
        D[i,j] = (-1)^(i+j) * cot((i-j)*T(π)/N) * jac_2
        D[i,i] -= D[i,j]
    end
    D
end




struct FourierSpectralViscosity{T<:Real, GridCompute, GridEvaluate, RFFT, IRFFT} <: AbstractDerivativeOperator{T}
    strength::T
    cutoff::Int
    D::FourierDerivativeOperator{T,GridCompute,GridEvaluate,RFFT,IRFFT}
end

function spectral_viscosity_operator(D::FourierDerivativeOperator{T},
                                     strength=T(1)/size(D,2),
                                     cutoff=round(Int, sqrt(size(D,2)))) where {T}
    FourierSpectralViscosity(strength, cutoff, D)
end

function Base.show(io::IO, Di::FourierSpectralViscosity{T}) where {T}
    grid = Di.D.grid_evaluate
    print(io, "Spectral viscosity operator for the periodic 1st derivative Fourier operator\n")
    print(io, "{T=", T, "} on a grid in [", first(grid), ", ", last(grid),
                "] using ", length(grid)-1, " modes\n")
    print(io, "with strength ε = ", Di.strength, " and cutoff m = ", Di.cutoff, ".\n")
end

grid(Di::FourierSpectralViscosity) = grid(Di.D)

function Base.A_mul_B!(dest::AbstractVector{T}, Di::FourierSpectralViscosity{T},
                        u::AbstractVector{T}) where {T}
    @unpack strength, cutoff, D = Di
    @unpack jac, tmp, rfft_plan, irfft_plan = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    A_mul_B!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(cutoff-1)
        tmp[j] = 0
    end
    @inbounds @simd for j in cutoff:(length(tmp))
        tmp[j] *= -strength * (j-1)^2 * jac^2 * exp(-((N-j+1)/(j-1-cutoff))^2)
    end
    #@inbounds tmp[end] = 0
    A_mul_B!(dest, irfft_plan, tmp)

    nothing
end
