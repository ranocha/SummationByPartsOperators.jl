
"""
    FourierDerivativeOperator{T<:Real, GridCompute, GridEvaluate, RFFT, BRFFT}

A derivative operator on a periodic grid with scalar type `T` computing the
first derivative using a spectral Fourier expansion via real discrete Fourier
transforms.
"""
struct FourierDerivativeOperator{T<:Real, Grid, RFFT, BRFFT} <: AbstractPeriodicDerivativeOperator{T}
    jac::T
    Δx::T
    grid_compute::Grid   # N-1 nodes, including the left and excluding the right boundary
    grid_evaluate::Grid #  N  nodes, including both boundaries
    tmp::Vector{Complex{T}}
    rfft_plan::RFFT
    brfft_plan::BRFFT

    function FourierDerivativeOperator(jac::T, Δx::T, grid_compute::Grid, grid_evaluate::Grid,
                                        tmp::Vector{Complex{T}}, rfft_plan::RFFT, brfft_plan::BRFFT) where {T<:Real, Grid, RFFT, BRFFT}
        @argcheck length(brfft_plan) == length(tmp) DimensionMismatch
        @argcheck length(brfft_plan) == (length(rfft_plan)÷2)+1 DimensionMismatch
        @argcheck length(grid_compute) == length(rfft_plan) DimensionMismatch
        @argcheck length(grid_compute) == length(grid_evaluate)-1 DimensionMismatch
        @argcheck first(grid_compute) == first(grid_evaluate)
        @argcheck step(grid_compute) ≈ step(grid_evaluate)
        @argcheck last(grid_compute) < last(grid_evaluate)

        new{T, Grid, RFFT, BRFFT}(jac, Δx, grid_compute, grid_evaluate, tmp, rfft_plan, brfft_plan)
    end
end

"""
    FourierDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}

Construct the `FourierDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` Fourier modes.
"""
function FourierDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}
    @argcheck N >= 1

    jac = 2*T(π) / (xmax - xmin) / N # / N because of brfft instead of BRFFT
    Δx = (xmax - xmin) / N
    grid_evaluate = linspace(xmin, xmax, N+1)
    grid_compute = linspace(xmin, grid_evaluate[end-1], N)
    u = zero.(grid_compute)
    rfft_plan = plan_rfft(u)
    uhat = rfft_plan*u
    brfft_plan = plan_brfft(uhat, N)

    FourierDerivativeOperator(jac, Δx, grid_compute, grid_evaluate, uhat, rfft_plan, brfft_plan)
end

function fourier_derivative_operator(xmin::T, xmax::T, N::Int) where {T<:Real}
    FourierDerivativeOperator(xmin, xmax, N)
end

derivative_order(D::FourierDerivativeOperator) = 1
Base.issymmetric(D::FourierDerivativeOperator) = false

function Base.show(io::IO, D::FourierDerivativeOperator{T}) where {T}
    grid = D.grid_evaluate
    print(io, "Periodic 1st derivative Fourier operator {T=", T, "} \n")
    print(io, "on a grid in [", first(grid), ", ", last(grid),
                "] using ", length(grid)-1, " modes. \n")
end


function Base.A_mul_B!(dest::AbstractVector{T}, D::FourierDerivativeOperator,
                        u::AbstractVector{T}) where {T}
    @unpack jac, tmp, rfft_plan, brfft_plan = D
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
    A_mul_B!(dest, brfft_plan, tmp)

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



"""
    FourierSpectralViscosity

A spectral viscosity operator on a periodic grid computing the derivative using
a spectral Fourier expansion via real discrete Fourier transforms.
"""
struct FourierSpectralViscosity{T<:Real, Grid, RFFT, BRFFT, SourceOfCoefficients} <: AbstractDerivativeOperator{T}
    strength::T
    cutoff::Int
    coefficients::Vector{T}
    D::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}
    source_of_coefficients::SourceOfCoefficients

    function FourierSpectralViscosity(source_of_coefficients::SourceOfCoefficients, 
                                        strength::T, cutoff::Int,
                                        D::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}) where {T<:Real, Grid, RFFT, BRFFT, SourceOfCoefficients}
        # precompute coefficients
        N = size(D, 1)
        jac = N * D.jac^2 # ^2: 2nd derivative; # *N: brfft instead of irfft
        cutoff = max(1, cutoff)
        coefficients = Array{T}(length(D.brfft_plan))
        set_filter_coefficients!(coefficients, source_of_coefficients, N, jac, strength, cutoff)
        new{T, Grid, RFFT, BRFFT, SourceOfCoefficients}(strength, cutoff, coefficients, D, source_of_coefficients)
    end
end

function spectral_viscosity_operator(source,
                                     D::FourierDerivativeOperator{T},
                                     strength=T(1)/size(D,2),
                                     cutoff=round(Int, sqrt(size(D,2)))) where {T}
    FourierSpectralViscosity(source, strength, cutoff, D)
end

@inline source_of_coefficients(Di::FourierSpectralViscosity) = (Di.source_of_coefficients)

function Base.show(io::IO, Di::FourierSpectralViscosity{T}) where {T}
    grid = Di.D.grid_evaluate
    print(io, "Spectral viscosity operator for the periodic 1st derivative Fourier operator\n")
    print(io, "{T=", T, "} on a grid in [", first(grid), ", ", last(grid),
                "] using ", length(grid)-1, " modes\n")
    print(io, "with strength ε = ", Di.strength, ", cutoff m = ", Di.cutoff, ", and coefficients from\n")
    print(io, "  ", Di.source_of_coefficients)
end

Base.issymmetric(Di::FourierSpectralViscosity) = true
grid(Di::FourierSpectralViscosity) = grid(Di.D)

function Base.A_mul_B!(dest::AbstractVector{T}, Di::FourierSpectralViscosity{T},
                        u::AbstractVector{T}) where {T}
    @unpack strength, cutoff, coefficients, D = Di
    @unpack jac, tmp, rfft_plan, brfft_plan = D
    N = size(D, 1)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
        @argcheck length(tmp) == length(coefficients)
    end

    A_mul_B!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp))
        tmp[j] *= coefficients[j]
    end
    A_mul_B!(dest, brfft_plan, tmp)

    nothing
end



"""
    Tadmor1989

Coefficients of the Fourier spectral viscosity given in
  Tadmor (1989)
  Convergence of Spectral Methods for Nonlinear Conservation Laws.
  SIAM Journal of Numerical Analysis 26, pp. 30-44.
"""
struct Tadmor1989 <: SourceOfCoefficients end

function Base.show(io::IO, ::Tadmor1989)
    print(io,
        "  Tadmor (1989) \n",
        "  Convergence of Spectral Methods for Nonlinear Conservation Laws. \n",
        "  SIAM Journal of Numerical Analysis 26, pp. 30-44. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                    source::Tadmor1989,
                                    N::Int, jac::T, strength::T, cutoff::Int) where {T<:Real}
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac * exp(-((N-j+1)/(j-1-cutoff))^2)
    end
end


"""
    Schochet1990

Coefficients of the Fourier spectral viscosity given in
  Schochet (1990)
  The Rate of Convergence of Spectral-Viscosity Methods for Periodic Scalar
    Conservation Laws.
  SIAM Journal of Numerical Analysis 27, pp. 1142-1159.
"""
struct Schochet1990 <: SourceOfCoefficients end

function Base.show(io::IO, ::Schochet1990)
    print(io,
        "  Schochet (1990) \n",
        "  The Rate of Convergence of Spectral-Viscosity Methods for Periodic Scalar\n",
        "    Conservation Laws. \n",
        "  SIAM Journal of Numerical Analysis 27, pp. 1142-1159. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                    source::Schochet1990,
                                    N::Int, jac::T, strength::T, cutoff::Int) where {T<:Real}
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff:min(2cutoff,length(coefficients))
        coefficients[j] = -strength * (j-1)^2 * jac * exp(-((2cutoff-j+1)/(j-1-cutoff))^2)
    end
    @inbounds @simd for j in 2cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac
    end
end
