
"""
    FourierDerivativeOperator2D{T<:Real}

A derivative operator on a two-dimensional periodic grid with scalar type `T`
computing the first derivatives using a spectral Fourier expansion via
real discrete Fourier transforms.
"""
struct FourierDerivativeOperator2D{T<:Real, Grid, RFFTx, BRFFTx, RFFTy, BRFFTy} <: AbstractPeriodicDerivativeOperator{T}
    jac_x::T
    jac_y::T
    Δx::T
    Δy::T
    grid_compute::Grid  # (Nx-1)x(Ny-1) nodes, including the left and excluding the right boundary
    grid_evaluate::Grid # (Nx)  x(Ny  ) nodes, including all boundaries
    tmp_x::Array{Complex{T}, 2}
    rfft_plan_x::RFFTx
    brfft_plan_x::BRFFTx
    tmp_y::Array{Complex{T}, 2}
    rfft_plan_y::RFFTy
    brfft_plan_y::BRFFTy

    function FourierDerivativeOperator2D(jac_x::T, jac_y::T, Δx::T, Δy::T, grid_compute::Grid, grid_evaluate::Grid,
                                         tmp_x::Array{Complex{T}, 2}, rfft_plan_x::RFFTx, brfft_plan_x::BRFFTx,
                                         tmp_y::Array{Complex{T}, 2}, rfft_plan_y::RFFTy, brfft_plan_y::BRFFTy) where {T<:Real, Grid, RFFTx, BRFFTx, RFFTy, BRFFTy}
        @argcheck size(brfft_plan_x) == size(tmp_x) DimensionMismatch
        @argcheck size(brfft_plan_x, 1) == (size(rfft_plan_x, 1)÷2)+1 DimensionMismatch
        @argcheck size(grid_compute) == size(rfft_plan_x) DimensionMismatch
        @argcheck size(grid_compute) == size(grid_evaluate).-1 DimensionMismatch

        new{T, Grid, RFFTx, BRFFTx, RFFTy, BRFFTy}(jac_x, jac_y, Δx, Δy, grid_compute, grid_evaluate,
                                                   tmp_x, rfft_plan_x, brfft_plan_x,
                                                   tmp_y, rfft_plan_y, brfft_plan_y)
    end
end

"""
    FourierDerivativeOperator2D(xmin, xmax, Nx, ymin, ymax, Ny)

Construct the `FourierDerivativeOperator` on a uniform grid between
`xmin` and `xmax` using `Nx` nodes and `ymin` and `ymax` using `Ny` nodes.
"""
function FourierDerivativeOperator2D(xmin::T, xmax::T, Nx::Int, ymin::T, ymax::T, Ny::Int) where {T<:Real}
    @argcheck Nx >= 1
    @argcheck Ny >= 1

    jac_x = 2*T(π) / (xmax - xmin) / Nx # / N because of brfft instead of BRFFT
    Δx = (xmax - xmin) / Nx
    jac_y = 2*T(π) / (ymax - ymin) / Ny # / N because of brfft instead of BRFFT
    Δy = (ymax - ymin) / Ny
    grid_evaluate = zeros(SVector{2,T}, Nx+1, Ny+1)
    grid_x = range(xmin, stop=xmax, length=Nx+1)
    grid_y = range(ymin, stop=ymax, length=Ny+1) # two boundary nodes
    for j in Base.OneTo(Ny+1), i in Base.OneTo(Nx+1)
        grid_evaluate[i,j] = SVector(grid_x[i], grid_y[j])
    end
    grid_compute  = grid_evaluate[1:end-1, 1:end-1]
    u = zero.(first.(grid_compute))
    rfft_plan_x = plan_rfft(u, 1)
    uhat_x = rfft_plan_x*u
    brfft_plan_x = plan_brfft(uhat_x, Nx, 1)
    rfft_plan_y = plan_rfft(u, 2)
    uhat_y = rfft_plan_y*u
    brfft_plan_y = plan_brfft(uhat_y, Ny, 2)

    FourierDerivativeOperator2D(jac_x, jac_y, Δx, Δy, grid_compute, grid_evaluate,
                                uhat_x, rfft_plan_x, brfft_plan_x,
                                uhat_y, rfft_plan_y, brfft_plan_y)
end

function fourier_derivative_operator(xmin::Real, xmax::Real, Nx::Int, ymin::Real, ymax::Real, Ny::Int)
    xmin, xmax, ymin, ymax = promote(xmin, xmax, ymin, ymax)
    FourierDerivativeOperator2D(xmin, xmax, Nx, ymin, ymax, Ny)
end

derivative_order(D::FourierDerivativeOperator2D) = 1
LinearAlgebra.issymmetric(D::FourierDerivativeOperator2D) = false

function Base.show(io::IO, D::FourierDerivativeOperator2D{T}) where {T}
    grid = D.grid_evaluate
    print(io, "Periodic 1st derivative Fourier operator {T=", T, "} in two space dimensions")
end


function mul!(dest::AbstractArray{T, 2}, D::FourierDerivativeOperator2D, u::AbstractArray{T, 2}, ::Val{:x}) where {T}
    @unpack jac_x, tmp_x, rfft_plan_x, brfft_plan_x = D
    Nx, Ny = size(D.grid_compute)
    @boundscheck begin
        @argcheck (Nx, Ny) == size(u)
        @argcheck (Nx, Ny) == size(dest)
    end

    mul!(tmp_x, rfft_plan_x, u)
    size_x, size_y = size(tmp_x)
    @inbounds @simd for j in Base.OneTo(size_y); for i in Base.OneTo(size_x-1)
        tmp_x[i, j] *= (i-1)*im * jac_x
    end end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    if iseven(Nx)
        @inbounds @simd for j in Base.OneTo(size_y)
            tmp_x[end, j] = zero(eltype(tmp_x))
        end
    else
        @inbounds @simd for j in Base.OneTo(size_y)
            tmp_x[end, j] *= (size_x-1)*im * jac_x
        end
    end
    mul!(dest, brfft_plan_x, tmp_x)
end

function mul!(dest::AbstractArray{T, 2}, D::FourierDerivativeOperator2D, u::AbstractArray{T, 2}, ::Val{:y}) where {T}
    @unpack jac_y, tmp_y, rfft_plan_y, brfft_plan_y = D
    Nx, Ny = size(D.grid_compute)
    @boundscheck begin
        @argcheck (Nx, Ny) == size(u)
        @argcheck (Nx, Ny) == size(dest)
    end

    mul!(tmp_y, rfft_plan_y, u)
    size_x, size_y = size(tmp_y)
    @inbounds @simd for j in Base.OneTo(size_y-1); for i in Base.OneTo(size_x)
        tmp_y[i, j] *= (j-1)*im * jac_y
    end end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    if iseven(Ny)
        @inbounds @simd for i in Base.OneTo(size_x)
            tmp_y[i, end] = zero(eltype(tmp_y))
        end
    else
        @inbounds @simd for i in Base.OneTo(size_x)
            tmp_y[i, end] *= (size_y-1)*im * jac_y
        end
    end
    mul!(dest, brfft_plan_y, tmp_y)
end

# TODO there is no 5 argument mul! in FFTW.jl...

function integrate(func, u::AbstractMatrix, D::FourierDerivativeOperator2D)
    @boundscheck begin
        length(u) == length(grid(D))
    end
    @unpack Δx, Δy = D

    @inbounds res = sum(func, u)

    (Δx * Δy) * res
end

function mass_matrix(D::FourierDerivativeOperator2D)
    @unpack Δx, Δy = D

    (Δx * Δy) * I
end
