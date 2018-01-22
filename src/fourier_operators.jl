
"""
    FourierDerivativeOperator{T<:Real, Grid, RFFT, BRFFT}

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
`xmax` using `N` nodes and `N÷2+1` complex Fourier modes.
"""
function FourierDerivativeOperator(xmin::T, xmax::T, N::Int) where {T<:Real}
    @argcheck N >= 1

    jac = 2*T(π) / (xmax - xmin) / N # / N because of brfft instead of BRFFT
    Δx = (xmax - xmin) / N
    grid_evaluate = linspace(xmin, xmax, N+1) # two boundary nodes
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
                "] using ", length(D.rfft_plan), " nodes and ",
                length(D.brfft_plan), " modes. \n")
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
basis using `N` nodes, see
Kopriva (2009) Implementing Spectral Methods for PDEs, Algorithm 18.
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



abstract type AbstractFourierViscosity{T} <: AbstractDerivativeOperator{T} end

@inline source_of_coefficients(Di::AbstractFourierViscosity) = (Di.source_of_coefficients)

Base.issymmetric(Di::AbstractFourierViscosity) = true
grid(Di::AbstractFourierViscosity) = grid(Di.D)

function Base.A_mul_B!(dest::AbstractVector{T}, Di::AbstractFourierViscosity{T},
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
    FourierSpectralViscosity

A spectral viscosity operator on a periodic grid computing the derivative using
a spectral Fourier expansion via real discrete Fourier transforms.
"""
struct FourierSpectralViscosity{T<:Real, Grid, RFFT, BRFFT, SourceOfCoefficients} <: AbstractFourierViscosity{T}
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
        jac = (N * D.jac)^2 / N # ^2: 2nd derivative; # *N / N: brfft instead of irfft
        cutoff = max(1, cutoff+1) # +1 because of 1 based indexing in Julia
        coefficients = Array{T}(length(D.brfft_plan))
        set_filter_coefficients!(coefficients, source_of_coefficients, jac, strength, cutoff)
        new{T, Grid, RFFT, BRFFT, SourceOfCoefficients}(strength, cutoff, coefficients, D, source_of_coefficients)
    end
end

"""
    spectral_viscosity_operator(source, D::FourierDerivativeOperator, 
                                strength=1/size(D,2), cutoff=sqrt(size(D,2)))

Construct the spectral viscosity operator for the Fourier derivative operator
`D` given in `source` with parameters `strength`, `cutoff`.
"""
function spectral_viscosity_operator(source,
                                     D::FourierDerivativeOperator{T},
                                     strength=T(1)/size(D,2),
                                     cutoff=round(Int, sqrt(size(D,2)))) where {T}
    FourierSpectralViscosity(source, strength, cutoff, D)
end

function Base.show(io::IO, Di::FourierSpectralViscosity{T}) where {T}
    grid = Di.D.grid_evaluate
    print(io, "Spectral viscosity operator for the periodic 1st derivative Fourier operator\n")
    print(io, "{T=", T, "} on a grid in [", first(grid), ", ", last(grid),
                "] using ", length(Di.D.rfft_plan), " nodes and ",
                length(Di.D.brfft_plan), " modes\n")
    print(io, "with strength ε = ", Di.strength, ", cutoff m = ", Di.cutoff, ", and coefficients from\n")
    print(io, Di.source_of_coefficients)
end



"""
    Tadmor1989

Coefficients of the Fourier spectral viscosity given in
  Tadmor (1989)
  Convergence of Spectral Methods for Nonlinear Conservation Laws.
  SIAM Journal on Numerical Analysis 26.1, pp. 30-44.
"""
struct Tadmor1989 <: SourceOfCoefficients end

function Base.show(io::IO, ::Tadmor1989)
    print(io,
        "  Tadmor (1989) \n",
        "  Convergence of Spectral Methods for Nonlinear Conservation Laws. \n",
        "  SIAM Journal on Numerical Analysis 26.1, pp. 30-44. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                    source::Tadmor1989,
                                    jac::T, strength::T, cutoff::Int) where {T<:Real}
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac
    end
end


"""
    MadayTadmor1989

Coefficients of the Fourier spectral viscosity given in
  Maday, Tadmor (1989)
  Analysis of the Spectral Vanishing Viscosity Method for Periodic Conservation
    Laws.
  SIAM Journal on Numerical Analysis 26.4, pp. 854-870.
"""
struct MadayTadmor1989 <: SourceOfCoefficients end

function Base.show(io::IO, ::MadayTadmor1989)
    print(io,
        "  Maday, Tadmor (1989) \n",
        "  Analysis of the Spectral Vanishing Viscosity Method for Periodic Conservation\n",
        "    Laws. \n",
        "  SIAM Journal on Numerical Analysis 26.4, pp. 854-870. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                    source::MadayTadmor1989,
                                    jac::T, strength::T, cutoff::Int) where {T<:Real}
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff:min(2cutoff,length(coefficients))
        coefficients[j] = -strength * (j-1)^2 * jac * (j-cutoff)/cutoff
    end
    @inbounds @simd for j in 2cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac
    end
end


"""
    TadmorWaagan2012Standard

Coefficients of the Fourier spectral viscosity given in
  Tadmor, Waagan (2012)
  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws.
  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009.
"""
struct TadmorWaagan2012Standard <: SourceOfCoefficients end

function Base.show(io::IO, ::TadmorWaagan2012Standard)
    print(io,
        "  Tadmor, Waagan (2012) \n",
        "  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws. \n",
        "  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                    source::TadmorWaagan2012Standard,
                                    jac::T, strength::T, cutoff::Int) where {T<:Real}
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff+1:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac * exp(-((length(coefficients)-j)/(j-cutoff))^2)
    end
end


"""
    TadmorWaagan2012Convergent

Coefficients of the Fourier spectral viscosity given in
  Tadmor, Waagan (2012)
  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws.
  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009.
See also
  Schochet (1990)
  The Rate of Convergence of Spectral-Viscosity Methods for Periodic Scalar
    Conservation Laws.
  SIAM Journal on Numerical Analysis 27.5, pp. 1142-1159.
"""
struct TadmorWaagan2012Convergent <: SourceOfCoefficients end

function Base.show(io::IO, ::TadmorWaagan2012Convergent)
    print(io,
        "  Tadmor, Waagan (2012) \n",
        "  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws. \n",
        "  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                  source::TadmorWaagan2012Convergent,
                                  jac::T, strength::T, cutoff::Int) where {T<:Real}
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff:min(2cutoff,length(coefficients))
        coefficients[j] = -strength * (j-1)^2 * jac * exp(-((2cutoff-j)/(j-cutoff))^2)
    end
    @inbounds @simd for j in 2cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac
    end
end


"""
    FourierSuperSpectralViscosity

A super spectral viscosity operator on a periodic grid computing the derivative
using a spectral Fourier expansion via real discrete Fourier transforms.
"""
struct FourierSuperSpectralViscosity{T<:Real, Grid, RFFT, BRFFT, SourceOfCoefficients} <: AbstractFourierViscosity{T}
    strength::T
    cutoff::Int
    order::Int
    coefficients::Vector{T}
    D::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}
    source_of_coefficients::SourceOfCoefficients

    function FourierSuperSpectralViscosity(source_of_coefficients::SourceOfCoefficients, 
                                        strength::T, cutoff::Int, order::Int,
                                        D::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}) where {T<:Real, Grid, RFFT, BRFFT, SourceOfCoefficients}
        # precompute coefficients
        @argcheck order >= 1

        N = size(D, 1)
        jac = (N * D.jac)^2order / N # ^2order: order-th derivative; # *N / N: brfft instead of irfft
        cutoff = max(1, cutoff+1) # +1 because of 1 based indexing in Julia
        coefficients = Array{T}(length(D.brfft_plan))
        set_filter_coefficients!(coefficients, source_of_coefficients, jac, strength, cutoff, order)
        new{T, Grid, RFFT, BRFFT, SourceOfCoefficients}(strength, cutoff, order, coefficients, D, source_of_coefficients)
    end
end

"""
    super_spectral_viscosity_operator(source, D::FourierDerivativeOperator, 
                                      order=1, 
                                      strength=1/size(D,2), 
                                      cutoff=sqrt(size(D,2)))

Construct the super spectral viscosity operator for the Fourier derivative
operator `D` given in `source` with parameters `order`, `strength`, `cutoff`.
"""
function super_spectral_viscosity_operator(source,
                                           D::FourierDerivativeOperator{T},
                                           order::Int=1,
                                           strength=T(1)/size(D,2)^(2order-1),
                                           cutoff::Int=round(Int, size(D,2)^(1-1/2order))) where {T}
    FourierSuperSpectralViscosity(source, strength, cutoff, order, D)
end

function Base.show(io::IO, Di::FourierSuperSpectralViscosity{T}) where {T}
    grid = Di.D.grid_evaluate
    print(io, "Super spectral viscosity operator for the periodic 1st derivative Fourier\n")
    print(io, "operator {T=", T, "} on a grid in [", first(grid), ", ", last(grid),
                "] using ", length(Di.D.rfft_plan), " nodes and ",
                length(Di.D.brfft_plan), " complex modes\n")
    print(io, "with strength ε = ", Di.strength, ", cutoff m = ", Di.cutoff,
              "order s = ", Di.order, ", and coefficients from\n")
    print(io, Di.source_of_coefficients)
end


"""
    Tadmor1993

Coefficients of the Fourier super spectral viscosity given in
  Tadmor (1993)
  Super Viscosity and Spectral Approximations of Nonlinear Conservation Laws.
  Numerical Methods for Fluid Dynamics IV, pp. 69-82.
"""
struct Tadmor1993 <: SourceOfCoefficients end

function Base.show(io::IO, ::Tadmor1993)
    print(io,
        "  Tadmor (1993) \n",
        "  Super Viscosity and Spectral Approximations of Nonlinear Conservation Laws. \n",
        "  Numerical Methods for Fluid Dynamics IV, pp. 69-82. \n")
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                    source::Tadmor1993,
                                    jac::T, strength::T, cutoff::Int, order::Int) where {T<:Real}
    @argcheck cutoff >= 1
    @argcheck order >= 1
    fill!(coefficients, zero(T))
    @inbounds @simd for j in cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2order * jac * (1 - (cutoff/j)^2order)
    end
end
