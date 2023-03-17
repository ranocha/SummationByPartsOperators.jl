
"""
    FourierDerivativeOperator{T}

A derivative operator on a periodic grid with real scalar type `T` computing the
first derivative using a spectral Fourier expansion via real discrete Fourier
transforms.
"""
@auto_hash_equals struct FourierDerivativeOperator{T<:Real, Grid, RFFT, BRFFT} <: AbstractPeriodicDerivativeOperator{T}
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
    FourierDerivativeOperator(xmin::T, xmax::T, N::Integer) where {T<:Real}

Construct the `FourierDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` nodes and `N÷2+1` complex Fourier modes.
"""
function FourierDerivativeOperator(xmin::T, xmax::T, N::Integer) where {T<:Real}
    @argcheck N >= 1

    jac = 2*T(π) / (xmax - xmin) / N # / N because of brfft instead of BRFFT
    Δx = (xmax - xmin) / N
    grid_evaluate = range(xmin, stop=xmax, length=N+1) # two boundary nodes
    grid_compute = range(xmin, stop=grid_evaluate[end-1], length=N)
    u = zero.(grid_compute)
    rfft_plan = plan_rfft(u)
    uhat = rfft_plan*u
    brfft_plan = plan_brfft(uhat, N)

    FourierDerivativeOperator(jac, Δx, grid_compute, grid_evaluate, uhat, rfft_plan, brfft_plan)
end

"""
    fourier_derivative_operator(xmin::Real, xmax::Real, N::Integer)
    fourier_derivative_operator(; xmin::Real, xmax::Real, N::Integer)

Construct the `FourierDerivativeOperator` on a uniform grid between `xmin` and
`xmax` using `N` nodes and `N÷2+1` complex Fourier modes.
"""
function fourier_derivative_operator(xmin::Real, xmax::Real, N::Integer)
    FourierDerivativeOperator(promote(xmin, xmax)..., N)
end

function fourier_derivative_operator(; xmin::Real, xmax::Real, N::Integer)
    FourierDerivativeOperator(promote(xmin, xmax)..., N)
end

derivative_order(D::FourierDerivativeOperator) = 1
LinearAlgebra.issymmetric(D::FourierDerivativeOperator) = false

function Base.show(io::IO, D::FourierDerivativeOperator)
    if get(io, :compact, false)
        print(io, "Periodic 1st derivative Fourier operator")
    else
        grid = D.grid_evaluate
        print(io, "Periodic 1st derivative Fourier operator {T=", eltype(D), "} \n")
        print(io, "on a grid in [", first(grid), ", ", last(grid),
                    "] using ", length(D.rfft_plan), " nodes and ",
                    length(D.brfft_plan), " modes")
    end
end


function mul!(dest::AbstractVector{T}, D::FourierDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack jac, tmp, rfft_plan, brfft_plan = D
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        tmp[j] *= (j-1)*im * jac
    end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    if iseven(N)
        @inbounds tmp[end] = zero(eltype(tmp))
    else
        @inbounds tmp[end] *= (length(tmp)-1)*im * jac
    end
    mul!(dest, brfft_plan, tmp)
end

# TODO there is no 5 argument mul! in FFTW.jl...
# function mul!(dest::AbstractVector{T}, D::FourierDerivativeOperator, u::AbstractVector{T}, α, β) where {T}
#     @unpack jac, tmp, rfft_plan, brfft_plan = D
#     N, _ = size(D)
#     @boundscheck begin
#         @argcheck N == length(u)
#         @argcheck N == length(dest)
#     end

#     mul!(tmp, rfft_plan, u)
#     @inbounds @simd for j in Base.OneTo(length(tmp)-1)
#         tmp[j] *= (j-1)*im * jac
#     end
#     if iseven(N)
#         @inbounds tmp[end] = zero(eltype(tmp))
#     else
#         @inbounds tmp[end] *= (length(tmp)-1)*im * jac
#     end
#     mul!(dest, brfft_plan, tmp, α, β)
# end



function integrate(func, u::AbstractVector, D::FourierDerivativeOperator)
    @boundscheck begin
        length(u) == length(grid(D))
    end
    @unpack Δx = D

    @inbounds res = sum(func, u)

    Δx * res
end

function mass_matrix(D::FourierDerivativeOperator)
    @unpack Δx = D

    Δx * I
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
    D = Array{T}(undef, N, N)
    @inbounds for j in 1:N, i in 1:N
        j == i && continue
        D[i,j] = (-1)^(i+j) * cot((i-j)*T(π)/N) * jac_2
        D[i,i] -= D[i,j]
    end
    D
end


function _coef_end(coef, D1::FourierDerivativeOperator{T}) where {T}
    if isodd(size(D1, 1))
        return coef
    else
        return ntuple(i -> isodd(i) ? coef[i] : zero(T), length(coef))
    end
end
# fallback for periodic FD operators etc.
_coef_end(coef, D1) = coef

@auto_hash_equals struct FourierPolynomialDerivativeOperator{T<:Real, Grid, RFFT, BRFFT, N} <: AbstractPeriodicDerivativeOperator{T}
    D1::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}
    coef::NTuple{N,T}
    coef_end::NTuple{N,T}

    function FourierPolynomialDerivativeOperator(D1::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}, coef::NTuple{N,T}) where {T<:Real, Grid, RFFT, BRFFT, N}
        coef_end = _coef_end(coef, D1)
        new{T,Grid,RFFT,BRFFT,N}(D1, coef, coef_end)
    end
end

function FourierPolynomialDerivativeOperator(D1::FourierDerivativeOperator)
    T = eltype(D1)
    FourierPolynomialDerivativeOperator(D1, (zero(T), one(T)))
end

Base.size(poly::FourierPolynomialDerivativeOperator) = size(poly.D1)
function LinearAlgebra.issymmetric(poly::FourierPolynomialDerivativeOperator)
    @unpack coef = poly
    all(iszero, coef[idx] for idx in eachindex(coef) if iseven(idx))
end
grid(poly::FourierPolynomialDerivativeOperator) = grid(poly.D1)

function Base.show(io::IO, poly::FourierPolynomialDerivativeOperator)
    if get(io, :compact, false)
        print(io, "Fourier polynomial")
    else
        print(io, "Fourier polynomial with coefficients\n")
        print(io, poly.coef)
        print(io, "\nof the operator:\n")
        print(io, poly.D1)
    end
end


function Base.:*(D1::FourierDerivativeOperator, D2::FourierDerivativeOperator)
    T = eltype(D1)
    @argcheck T == eltype(D2) ArgumentError
    @argcheck D1.jac == D2.jac ArgumentError
    @argcheck D1.Δx == D2.Δx ArgumentError
    @argcheck D1.grid_compute == D2.grid_compute DimensionMismatch
    @argcheck D1.grid_evaluate == D2.grid_evaluate DimensionMismatch

    FourierPolynomialDerivativeOperator(D1, (zero(T), zero(T), one(T)))
end

function Base.literal_pow(::typeof(^), D1::FourierDerivativeOperator, ::Val{P}) where {P}
    T = eltype(D1)
    coef = Base.setindex( ntuple(_->zero(T), Val{P+1}()), one(T), P+1)
    FourierPolynomialDerivativeOperator(D1, coef)
end

function Base.:*(factor::Union{Real,Integer}, poly::FourierPolynomialDerivativeOperator)
    @unpack coef = poly
    for idx in 1:length(coef)
        coef = Base.setindex(coef, factor*coef[idx], idx)
    end

    FourierPolynomialDerivativeOperator(poly.D1, coef)
end

function Base.:*(poly::FourierPolynomialDerivativeOperator, factor::Union{Real,Integer})
    factor * poly
end

function Base.:*(D::FourierDerivativeOperator, factor::Union{Real,Integer})
    FourierPolynomialDerivativeOperator(D) * factor
end

function Base.:*(factor::Union{Real,Integer}, D::FourierDerivativeOperator)
    D * factor
end

function Base.:*(poly::FourierPolynomialDerivativeOperator, scaling::UniformScaling)
    scaling.λ * poly
end

function Base.:*(scaling::UniformScaling, poly::FourierPolynomialDerivativeOperator)
    poly * scaling
end

function Base.:*(D::FourierDerivativeOperator, scaling::UniformScaling)
    scaling * FourierPolynomialDerivativeOperator(D)
end

function Base.:*(scaling::UniformScaling, D::FourierDerivativeOperator)
    FourierPolynomialDerivativeOperator(D) * scaling
end

function Base.:*(poly1::FourierPolynomialDerivativeOperator, poly2::FourierPolynomialDerivativeOperator)
    T = eltype(poly1.D1)
    @argcheck T == eltype(poly2.D1) ArgumentError
    @argcheck poly1.D1.jac == poly2.D1.jac ArgumentError
    @argcheck poly1.D1.Δx == poly2.D1.Δx ArgumentError
    @argcheck poly1.D1.grid_compute == poly2.D1.grid_compute DimensionMismatch
    @argcheck poly1.D1.grid_evaluate == poly2.D1.grid_evaluate DimensionMismatch

    coef = mul_poly(poly1.coef, poly2.coef)
    FourierPolynomialDerivativeOperator(poly1.D1, coef)
end

function Base.:*(poly1::FourierPolynomialDerivativeOperator, D2::FourierDerivativeOperator)
    poly1 * FourierPolynomialDerivativeOperator(D2)
end

function Base.:*(D1::FourierDerivativeOperator, poly2::FourierPolynomialDerivativeOperator)
    FourierPolynomialDerivativeOperator(D1) * poly2
end


function Base.:+(poly1::FourierPolynomialDerivativeOperator, poly2::FourierPolynomialDerivativeOperator)
    T = eltype(poly1.D1)
    @argcheck T == eltype(poly2.D1) ArgumentError
    @argcheck poly1.D1.jac == poly2.D1.jac ArgumentError
    @argcheck poly1.D1.Δx == poly2.D1.Δx ArgumentError
    @argcheck poly1.D1.grid_compute == poly2.D1.grid_compute DimensionMismatch
    @argcheck poly1.D1.grid_evaluate == poly2.D1.grid_evaluate DimensionMismatch

    coef = add_poly(poly1.coef, poly2.coef)
    FourierPolynomialDerivativeOperator(poly1.D1, coef)
end

function Base.:+(D1::FourierDerivativeOperator, poly2::FourierPolynomialDerivativeOperator)
    FourierPolynomialDerivativeOperator(D1) + poly2
end

function Base.:+(poly1::FourierPolynomialDerivativeOperator, D2::FourierDerivativeOperator)
    poly1 + FourierPolynomialDerivativeOperator(D2)
end

function Base.:+(poly::FourierPolynomialDerivativeOperator, scaling::UniformScaling)
    @unpack coef = poly
    coef = Base.setindex(coef, coef[1] + scaling.λ, 1)

    FourierPolynomialDerivativeOperator(poly.D1, coef)
end

function Base.:+(scaling::UniformScaling, poly::FourierPolynomialDerivativeOperator)
    poly + scaling
end

function Base.:+(D::FourierDerivativeOperator, scaling::UniformScaling)
    FourierPolynomialDerivativeOperator(D) + scaling
end

function Base.:+(scaling::UniformScaling, D::FourierDerivativeOperator)
    D + scaling
end


function Base.:-(poly1::FourierPolynomialDerivativeOperator, poly2::FourierPolynomialDerivativeOperator)
    T = eltype(poly1.D1)
    @argcheck T == eltype(poly2.D1) ArgumentError
    @argcheck poly1.D1.jac == poly2.D1.jac ArgumentError
    @argcheck poly1.D1.Δx == poly2.D1.Δx ArgumentError
    @argcheck poly1.D1.grid_compute == poly2.D1.grid_compute DimensionMismatch
    @argcheck poly1.D1.grid_evaluate == poly2.D1.grid_evaluate DimensionMismatch

    coef = subtract_poly(poly1.coef, poly2.coef)
    FourierPolynomialDerivativeOperator(poly1.D1, coef)
end

function Base.:-(D1::FourierDerivativeOperator, poly2::FourierPolynomialDerivativeOperator)
    FourierPolynomialDerivativeOperator(D1) - poly2
end

function Base.:-(poly1::FourierPolynomialDerivativeOperator, D2::FourierDerivativeOperator)
    poly1 - FourierPolynomialDerivativeOperator(D2)
end

function Base.:-(poly::FourierPolynomialDerivativeOperator, scaling::UniformScaling)
    @unpack coef = poly
    coef = Base.setindex(coef, coef[1] - scaling.λ, 1)

    FourierPolynomialDerivativeOperator(poly.D1, coef)
end

function Base.:-(scaling::UniformScaling, poly::FourierPolynomialDerivativeOperator)
    @unpack coef = poly
    coef = Base.setindex(coef, scaling.λ - coef[1], 1)
    for idx in 2:length(coef)
        coef = Base.setindex(coef, -coef[idx], idx)
    end

    FourierPolynomialDerivativeOperator(poly.D1, coef)
end

function Base.:-(D::FourierDerivativeOperator, scaling::UniformScaling)
    FourierPolynomialDerivativeOperator(D) - scaling
end

function Base.:-(scaling::UniformScaling, D::FourierDerivativeOperator)
    scaling - D
end


function mul!(dest::AbstractVector{T}, poly::FourierPolynomialDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack D1, coef, coef_end = poly
    @unpack jac, tmp, rfft_plan, brfft_plan = D1
    N, _ = size(D1)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        # *N ) / N: brfft instead of irfft
        tmp[j] *= evalpoly((j-1)*im * jac*N, coef) / N
    end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    @inbounds tmp[end] *= evalpoly((length(tmp)-1)*im * jac*N, coef_end) / N
    mul!(dest, brfft_plan, tmp)
end

function LinearAlgebra.ldiv!(dest::AbstractVector{T}, rat::FourierPolynomialDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack D1, coef, coef_end = rat
    @unpack jac, tmp, rfft_plan, brfft_plan = D1
    N, _ = size(D1)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        # *N ) / N: brfft instead of irfft
        tmp[j] /= (evalpoly((j-1)*im * jac*N, coef) * N)
    end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    @inbounds tmp[end] /= (evalpoly((length(tmp)-1)*im * jac*N, coef_end) * N)
    mul!(dest, brfft_plan, tmp)
end



@auto_hash_equals struct FourierRationalDerivativeOperator{T<:Real, Grid, RFFT, BRFFT, Nnum, Nden} <: AbstractPeriodicDerivativeOperator{T}
    D1::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}
    num_coef::NTuple{Nnum,T}
    num_coef_end::NTuple{Nnum,T}
    den_coef::NTuple{Nden,T}
    den_coef_end::NTuple{Nden,T}

    function FourierRationalDerivativeOperator(D1::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}, num_coef::NTuple{Nnum,T}, den_coef::NTuple{Nden,T}) where {T<:Real, Grid, RFFT, BRFFT, Nnum, Nden}
        num_coef_end = _coef_end(num_coef, D1)
        den_coef_end = _coef_end(den_coef, D1)
        new{T,Grid,RFFT,BRFFT,Nnum,Nden}(D1, num_coef, num_coef_end, den_coef, den_coef_end)
    end
end

Base.size(rat::FourierRationalDerivativeOperator) = size(rat.D1)
grid(rat::FourierRationalDerivativeOperator) = grid(rat.D1)

function Base.show(io::IO, rat::FourierRationalDerivativeOperator)
    if get(io, :compact, false)
        print(io, "Rational Fourier operator")
    else
        print(io, "Rational Fourier operator with coefficients\n")
        print(io, rat.num_coef)
        print(io, "\nand\n")
        print(io, rat.den_coef)
        print(io, "\nof the operator:\n")
        print(io, rat.D1)
    end
end

function LinearAlgebra.issymmetric(rat::FourierRationalDerivativeOperator)
    @unpack num_coef, den_coef = rat
    num_is_even = all(iszero, num_coef[idx] for idx in eachindex(num_coef) if iseven(idx))
    den_is_even = all(iszero, den_coef[idx] for idx in eachindex(den_coef) if iseven(idx))

    num_is_even == den_is_even
end


function FourierRationalDerivativeOperator(num::FourierPolynomialDerivativeOperator)
    T = eltype(num)
    FourierRationalDerivativeOperator(num.D1, num.coef, (one(T),))
end

function FourierRationalDerivativeOperator(D::FourierDerivativeOperator)
    FourierRationalDerivativeOperator(FourierPolynomialDerivativeOperator(D))
end

function Base.:/(num::FourierPolynomialDerivativeOperator, den::FourierPolynomialDerivativeOperator)
    @argcheck num.D1.jac == den.D1.jac ArgumentError
    @argcheck num.D1.Δx == den.D1.Δx ArgumentError
    @argcheck num.D1.grid_compute == den.D1.grid_compute DimensionMismatch
    @argcheck num.D1.grid_evaluate == den.D1.grid_evaluate DimensionMismatch

    FourierRationalDerivativeOperator(num.D1, num.coef, den.coef)
end

function Base.:/(num::FourierDerivativeOperator, den::FourierPolynomialDerivativeOperator)
    FourierPolynomialDerivativeOperator(num) / den
end

function Base.:/(num::FourierPolynomialDerivativeOperator, den::FourierDerivativeOperator)
    num / FourierPolynomialDerivativeOperator(den)
end

function Base.inv(rat::FourierRationalDerivativeOperator)
    FourierRationalDerivativeOperator(rat.D1, rat.den_coef, rat.num_coef)
end

function Base.inv(den::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator})
    inv(FourierRationalDerivativeOperator(den))
end


function Base.:+(rat1::FourierRationalDerivativeOperator, rat2::FourierRationalDerivativeOperator)
    T = eltype(rat1)
    @argcheck T == eltype(rat2) ArgumentError
    @argcheck rat1.D1.jac == rat2.D1.jac ArgumentError
    @argcheck rat1.D1.Δx == rat2.D1.Δx ArgumentError
    @argcheck rat1.D1.grid_compute == rat2.D1.grid_compute DimensionMismatch
    @argcheck rat1.D1.grid_evaluate == rat2.D1.grid_evaluate DimensionMismatch

    num_coef = add_poly(mul_poly(rat1.num_coef, rat2.den_coef), mul_poly(rat1.den_coef, rat2.num_coef))
    den_coef = mul_poly(rat1.den_coef, rat2.den_coef)
    FourierRationalDerivativeOperator(rat1.D1, num_coef, den_coef)
end

function Base.:+(rat1::FourierRationalDerivativeOperator, rat2::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator})
    rat1 + FourierRationalDerivativeOperator(rat2)
end

function Base.:+(rat1::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator}, rat2::FourierRationalDerivativeOperator)
    FourierRationalDerivativeOperator(rat1) + rat2
end

function Base.:-(rat1::FourierRationalDerivativeOperator, rat2::FourierRationalDerivativeOperator)
    T = eltype(rat1)
    @argcheck T == eltype(rat2) ArgumentError
    @argcheck rat1.D1.jac == rat2.D1.jac ArgumentError
    @argcheck rat1.D1.Δx == rat2.D1.Δx ArgumentError
    @argcheck rat1.D1.grid_compute == rat2.D1.grid_compute DimensionMismatch
    @argcheck rat1.D1.grid_evaluate == rat2.D1.grid_evaluate DimensionMismatch

    num_coef = subtract_poly(mul_poly(rat1.num_coef, rat2.den_coef), mul_poly(rat1.den_coef, rat2.num_coef))
    den_coef = mul_poly(rat1.den_coef, rat2.den_coef)
    FourierRationalDerivativeOperator(rat1.D1, num_coef, den_coef)
end

function Base.:-(rat1::FourierRationalDerivativeOperator, rat2::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator})
    rat1 - FourierRationalDerivativeOperator(rat2)
end

function Base.:-(rat1::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator}, rat2::FourierRationalDerivativeOperator)
    FourierRationalDerivativeOperator(rat1) - rat2
end

function Base.:*(rat1::FourierRationalDerivativeOperator, rat2::FourierRationalDerivativeOperator)
    T = eltype(rat1)
    @argcheck T == eltype(rat2) ArgumentError
    @argcheck rat1.D1.jac == rat2.D1.jac ArgumentError
    @argcheck rat1.D1.Δx == rat2.D1.Δx ArgumentError
    @argcheck rat1.D1.grid_compute == rat2.D1.grid_compute DimensionMismatch
    @argcheck rat1.D1.grid_evaluate == rat2.D1.grid_evaluate DimensionMismatch

    num_coef = mul_poly(rat1.num_coef, rat2.num_coef)
    den_coef = mul_poly(rat1.den_coef, rat2.den_coef)
    FourierRationalDerivativeOperator(rat1.D1, num_coef, den_coef)
end

function Base.:*(rat1::FourierRationalDerivativeOperator, rat2::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator})
    rat1 * FourierRationalDerivativeOperator(rat2)
end

function Base.:*(rat1::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator}, rat2::FourierRationalDerivativeOperator)
    FourierRationalDerivativeOperator(rat1) * rat2
end

function Base.:/(rat1::FourierRationalDerivativeOperator, rat2::FourierRationalDerivativeOperator)
    T = eltype(rat1)
    @argcheck T == eltype(rat2) ArgumentError
    @argcheck rat1.D1.jac == rat2.D1.jac ArgumentError
    @argcheck rat1.D1.Δx == rat2.D1.Δx ArgumentError
    @argcheck rat1.D1.grid_compute == rat2.D1.grid_compute DimensionMismatch
    @argcheck rat1.D1.grid_evaluate == rat2.D1.grid_evaluate DimensionMismatch

    num_coef = mul_poly(rat1.num_coef, rat2.den_coef)
    den_coef = mul_poly(rat1.den_coef, rat2.num_coef)
    FourierRationalDerivativeOperator(rat1.D1, num_coef, den_coef)
end

function Base.:/(rat1::FourierRationalDerivativeOperator, rat2::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator})
    rat1 / FourierRationalDerivativeOperator(rat2)
end

function Base.:/(rat1::Union{FourierDerivativeOperator,FourierPolynomialDerivativeOperator}, rat2::FourierRationalDerivativeOperator)
    FourierRationalDerivativeOperator(rat1) / rat2
end


function mul!(dest::AbstractVector{T}, rat::FourierRationalDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack D1, num_coef, num_coef_end, den_coef, den_coef_end = rat
    @unpack jac, tmp, rfft_plan, brfft_plan = D1
    N, _ = size(D1)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        # *N ) / N: brfft instead of irfft
        tmp[j] *= evalpoly((j-1)*im * jac*N, num_coef) / (N * evalpoly((j-1)*im * jac*N, den_coef))
    end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    @inbounds tmp[end] *= evalpoly((length(tmp)-1)*im * jac*N, num_coef_end) / (N * evalpoly((length(tmp)-1)*im * jac*N, den_coef_end))
    mul!(dest, brfft_plan, tmp)
end

function LinearAlgebra.ldiv!(dest::AbstractVector{T}, rat::FourierRationalDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack D1, num_coef, num_coef_end, den_coef, den_coef_end = rat
    @unpack jac, tmp, rfft_plan, brfft_plan = D1
    N, _ = size(D1)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        # *N ) / N: brfft instead of irfft
        tmp[j] *= evalpoly((j-1)*im * jac*N, den_coef) / (N * evalpoly((j-1)*im * jac*N, num_coef))
    end
    # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
    @inbounds tmp[end] *= evalpoly((N-1)*im * jac*N, den_coef_end) / (N * evalpoly((N-1)*im * jac*N, num_coef_end))
    mul!(dest, brfft_plan, tmp)
end

function Base.:\(rat::Union{FourierRationalDerivativeOperator,FourierPolynomialDerivativeOperator}, u::AbstractVector{T}) where {T}
    dest = similar(u)
    ldiv!(dest, rat, u)
end




@auto_hash_equals struct PeriodicDerivativeOperatorQuotient{T<:Real, numDtype<:Union{PeriodicDerivativeOperator{T},FourierDerivativeOperator{T}}, denDtype<:Union{PeriodicDerivativeOperator{T},FourierDerivativeOperator{T}}, Nnum, Nden, RFFT, IRFFT} <: AbstractPeriodicDerivativeOperator{T}
    num_D::numDtype
    den_D::denDtype
    num_coef::NTuple{Nnum,T}
    den_coef::NTuple{Nden,T}
    num_coef_end::NTuple{Nnum,T}
    den_coef_end::NTuple{Nden,T}
    tmp::Vector{Complex{T}}
    num_eigval::Vector{Complex{T}}
    den_eigval::Vector{Complex{T}}
    rfft_plan::RFFT
    irfft_plan::IRFFT

    function PeriodicDerivativeOperatorQuotient(num_D::numDtype, den_D::denDtype, num_coef::NTuple{Nnum,T}, den_coef::NTuple{Nden,T}, tmp::Vector{Complex{T}}, num_eigval::Vector{Complex{T}}, den_eigval::Vector{Complex{T}}, rfft_plan::RFFT, irfft_plan::IRFFT) where {T<:Real, numDtype<:Union{PeriodicDerivativeOperator{T},FourierDerivativeOperator{T}}, denDtype<:Union{PeriodicDerivativeOperator{T},FourierDerivativeOperator{T}}, Nnum, Nden, RFFT, IRFFT}
        @argcheck length(irfft_plan) == length(tmp) DimensionMismatch
        @argcheck length(irfft_plan) == length(num_eigval) DimensionMismatch
        @argcheck length(irfft_plan) == length(den_eigval) DimensionMismatch
        @argcheck length(irfft_plan) == (length(rfft_plan)÷2)+1 DimensionMismatch
        @argcheck length(grid(num_D)) == length(rfft_plan) DimensionMismatch
        @argcheck length(grid(den_D)) == length(rfft_plan) DimensionMismatch

        num_coef_end = _coef_end(num_coef, num_D)
        den_coef_end = _coef_end(den_coef, den_D)

        new{T,numDtype,denDtype,Nnum,Nden,RFFT,IRFFT}(num_D, den_D, num_coef, den_coef, num_coef_end, den_coef_end, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
    end
end

Base.size(quot::PeriodicDerivativeOperatorQuotient) = size(quot.num_D)
grid(quot::PeriodicDerivativeOperatorQuotient) = grid(quot.num_D)

function Base.show(io::IO, quot::PeriodicDerivativeOperatorQuotient)
    if get(io, :compact, false)
        print(io, "Rational operator (", quot.num_D, ", ", quot.den_D, ")")
    else
        print(io, "Quotient of the polynomial with coefficients\n")
        print(io, quot.num_coef)
        print(io, "\nof the operator:\n")
        print(io, quot.num_D)
        print(io, "and the polynomial with coefficients\n")
        print(io, quot.den_coef)
        print(io, "\nof the operator:\n")
        print(io, quot.den_D)
    end
end


function _eigvals!(eigval::Vector{Complex{T}}, D::FourierDerivativeOperator{T}) where {T<:Real}
    N = length(grid(D))
    jac = D.jac * N
    @inbounds for idx in 1:(length(eigval))
        eigval[idx] = (idx-1)*im * jac
    end

    eigval
end

function Base.://(num::FourierDerivativeOperator, den::PeriodicRationalDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError
    @argcheck den.den_coef == (one(T),) ArgumentError

    @unpack tmp, rfft_plan = num
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D)

    PeriodicDerivativeOperatorQuotient(num, den.D, (zero(T), one(T)), den.num_coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::FourierPolynomialDerivativeOperator, den::PeriodicRationalDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError
    @argcheck den.den_coef == (one(T),) ArgumentError

    @unpack tmp, rfft_plan = num.D1
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num.D1)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D)

    PeriodicDerivativeOperatorQuotient(num.D1, den.D, num.coef, den.num_coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::PeriodicDerivativeOperator, den::PeriodicRationalDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError
    @argcheck den.den_coef == (one(T),) ArgumentError

    x = grid(num)
    u = zero.(x)
    rfft_plan = plan_rfft(u)
    tmp = rfft_plan * u
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D)

    PeriodicDerivativeOperatorQuotient(num, den.D, (zero(T), one(T)), den.num_coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::PeriodicRationalDerivativeOperator, den::PeriodicRationalDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError
    @argcheck num.den_coef == (one(T),) ArgumentError
    @argcheck den.den_coef == (one(T),) ArgumentError

    x = grid(num)
    u = zero.(x)
    rfft_plan = plan_rfft(u)
    tmp = rfft_plan * u
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num.D)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D)

    PeriodicDerivativeOperatorQuotient(num.D, den.D, num.num_coef, den.num_coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::FourierDerivativeOperator, den::FourierPolynomialDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError

    @unpack tmp, rfft_plan = den.D1
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D1)

    PeriodicDerivativeOperatorQuotient(num, den.D1, (zero(T), one(T)), den.coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::FourierPolynomialDerivativeOperator, den::FourierPolynomialDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError

    @unpack tmp, rfft_plan = den.D1
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num.D1)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D1)

    PeriodicDerivativeOperatorQuotient(num.D1, den.D1, num.coef, den.coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::PeriodicDerivativeOperator, den::FourierPolynomialDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError

    @unpack tmp, rfft_plan = den.D1
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D1)

    PeriodicDerivativeOperatorQuotient(num, den.D1, (zero(T), one(T)), den.coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end

function Base.://(num::PeriodicRationalDerivativeOperator, den::FourierPolynomialDerivativeOperator)
    T = eltype(num)
    @argcheck T == eltype(den) ArgumentError
    @argcheck grid(num) == grid(den) ArgumentError
    @argcheck num.den_coef == (one(T),) ArgumentError

    @unpack tmp, rfft_plan = den.D1
    irfft_plan = plan_irfft(tmp, length(grid(num)))

    num_eigval = similar(tmp)
    _eigvals!(num_eigval, num.D)

    den_eigval = similar(tmp)
    _eigvals!(den_eigval, den.D1)

    PeriodicDerivativeOperatorQuotient(num.D, den.D1, num.num_coef, den.coef, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan)
end


function mul!(dest::AbstractVector{T}, quot::PeriodicDerivativeOperatorQuotient, u::AbstractVector{T}) where {T}
    @unpack num_D, num_coef, den_coef, num_coef_end, den_coef_end, tmp, num_eigval, den_eigval, rfft_plan, irfft_plan = quot
    N, _ = size(num_D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp)-1)
        tmp[j] *= evalpoly(num_eigval[j], num_coef) / evalpoly(den_eigval[j], den_coef)
    end
    tmp[end] *= evalpoly(num_eigval[end], num_coef_end) / evalpoly(den_eigval[end], den_coef_end)
    mul!(dest, irfft_plan, tmp)
end




"""
    ConstantFilter(D::FourierDerivativeOperator, filter)

Create a modal filter with constant parameters adapted to the Fourier
derivative operator `D` with parameters given by the filter function `filter`.
"""
function ConstantFilter(D::FourierDerivativeOperator{T}, filter) where {T}
    Np1 = length(D.brfft_plan)
    coefficients = Array{T}(undef, Np1)
    set_filter_coefficients!(coefficients, filter)
    tmp = copy(D.tmp)
    modal2nodal = plan_irfft(D.tmp, length(D.rfft_plan))
    nodal2modal = D.rfft_plan
    ConstantFilter(coefficients, nodal2modal, modal2nodal, tmp, filter)
end



abstract type AbstractFourierViscosity{T} <: AbstractPeriodicDerivativeOperator{T} end

@inline source_of_coefficients(Di::AbstractFourierViscosity) = (Di.source_of_coefficients)

LinearAlgebra.issymmetric(Di::AbstractFourierViscosity) = true
grid(Di::AbstractFourierViscosity) = grid(Di.D)

function mul!(dest::AbstractVector{T}, Di::AbstractFourierViscosity{T}, u::AbstractVector{T}) where {T}
    @unpack coefficients, D = Di
    @unpack tmp, rfft_plan, brfft_plan = D
    N = size(D, 1)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
        @argcheck length(tmp) == length(coefficients)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp))
        tmp[j] *= coefficients[j]
    end
    mul!(dest, brfft_plan, tmp)

    nothing
end


"""
    FourierConstantViscosity

Fourier viscosity operator with constant coefficients for the periodic 1st
derivative Fourier operator.
"""
@auto_hash_equals struct FourierConstantViscosity{T<:Real, Grid, RFFT, BRFFT} <: AbstractFourierViscosity{T}
    coefficients::Vector{T}
    D::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}
    parameters
    source_of_coefficients

    function FourierConstantViscosity(D::FourierDerivativeOperator{T,Grid,RFFT,BRFFT}, parameters, source_of_coefficients) where {T<:Real, Grid, RFFT, BRFFT}
        # precompute coefficients
        N = size(D,1)
        jac = N * D.jac # *N: brfft instead of irfft
        coefficients = Vector{T}(undef, length(D.brfft_plan))
        set_filter_coefficients!(coefficients, jac, N, parameters, source_of_coefficients)
        new{T,Grid,RFFT,BRFFT}(coefficients, D, parameters, source_of_coefficients)
    end
end

function Base.show(io::IO, Di::FourierConstantViscosity)
    if get(io, :compact, false)
        print(io, "Fourier viscosity operator with constant coefficients of ")
    else
        grid = Di.D.grid_evaluate
        print(io, "Fourier viscosity operator with constant coefficients for the periodic 1st\n")
        print(io, "derivative Fourier operator {T=", eltype(Di), "} on a grid in [",
                first(grid), ", ", last(grid), "]\n")
        print(io, "using ", length(Di.D.rfft_plan), " nodes and ",
                length(Di.D.brfft_plan), " modes with coefficients of ")
    end
    print(io, Di.source_of_coefficients)
end


function dissipation_operator(source_of_coefficients, D::FourierDerivativeOperator; kwargs...)
    parameters = get_parameters(source_of_coefficients, D; kwargs...)
    FourierConstantViscosity(D, parameters, source_of_coefficients)
end


"""
    Tadmor1989()

Coefficients of the Fourier spectral viscosity given in
- Tadmor (1989)
  Convergence of Spectral Methods for Nonlinear Conservation Laws.
  SIAM Journal on Numerical Analysis 26.1, pp. 30-44.
"""
struct Tadmor1989 <: SourceOfCoefficients end

function Base.show(io::IO, source::Tadmor1989)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
            "Tadmor (1989) \n",
            "  Convergence of Spectral Methods for Nonlinear Conservation Laws. \n",
            "  SIAM Journal on Numerical Analysis 26.1, pp. 30-44. \n")
    end
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                  jac::T, N::Int,
                                  parameters, source::Tadmor1989) where {T<:Real}
    @unpack strength, cutoff = parameters
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    jac = jac^2 / N # ^2: 2nd derivative; /N: brfft instead of irfft
    @simd for j in cutoff:length(coefficients)
        @inbounds coefficients[j] = -strength * (j-1)^2 * jac
    end
end


"""
    MadayTadmor1989()

Coefficients of the Fourier spectral viscosity given in
- Maday, Tadmor (1989)
  Analysis of the Spectral Vanishing Viscosity Method for Periodic Conservation
    Laws.
  SIAM Journal on Numerical Analysis 26.4, pp. 854-870.
"""
struct MadayTadmor1989 <: SourceOfCoefficients end

function Base.show(io::IO, source::MadayTadmor1989)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
            "Maday, Tadmor (1989) \n",
            "  Analysis of the Spectral Vanishing Viscosity Method for Periodic Conservation\n",
            "    Laws. \n",
            "  SIAM Journal on Numerical Analysis 26.4, pp. 854-870.")
    end
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                  jac::T, N::Int,
                                  parameters, source::MadayTadmor1989) where {T<:Real}
    @unpack strength, cutoff = parameters
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    jac = jac^2 / N # ^2: 2nd derivative; /N: brfft instead of irfft
    @inbounds @simd for j in cutoff:min(2cutoff,length(coefficients))
        coefficients[j] = -strength * (j-1)^2 * jac * (j-cutoff)/cutoff
    end
    @inbounds @simd for j in 2cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac
    end
end


"""
    TadmorWaagan2012Standard()

Coefficients of the Fourier spectral viscosity given in
- Tadmor, Waagan (2012)
  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws.
  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009.
"""
struct TadmorWaagan2012Standard <: SourceOfCoefficients end

function Base.show(io::IO, source::TadmorWaagan2012Standard)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
            "Tadmor, Waagan (2012) \n",
            "  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws. \n",
            "  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009.")
    end
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                  jac::T, N::Int,
                                  parameters, source::TadmorWaagan2012Standard) where {T<:Real}
    @unpack strength, cutoff = parameters
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    jac = jac^2 / N # ^2: 2nd derivative; /N: brfft instead of irfft
    @inbounds @simd for j in cutoff+1:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac * exp(-((length(coefficients)-j)/(j-cutoff))^2)
    end
end


"""
    TadmorWaagan2012Convergent()

Coefficients of the Fourier spectral viscosity given in
- Tadmor, Waagan (2012)
  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws.
  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009.

See also
- Schochet (1990)
  The Rate of Convergence of Spectral-Viscosity Methods for Periodic Scalar
    Conservation Laws.
  SIAM Journal on Numerical Analysis 27.5, pp. 1142-1159.
"""
struct TadmorWaagan2012Convergent <: SourceOfCoefficients end

function Base.show(io::IO, source::TadmorWaagan2012Convergent)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
            "Tadmor, Waagan (2012) \n",
            "  Adaptive Spectral Viscosity for Hyperbolic Conservation Laws. \n",
            "  SIAM Journal on Scientific Computing 34.2, pp. A993-A1009.")
    end
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                  jac::T, N::Int,
                                  parameters, source::TadmorWaagan2012Convergent) where {T<:Real}
    @unpack strength, cutoff = parameters
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    jac = jac^2 / N # ^2: 2nd derivative; /N: brfft instead of irfft
    @inbounds @simd for j in cutoff:min(2cutoff,length(coefficients))
        coefficients[j] = -strength * (j-1)^2 * jac * exp(-((2cutoff-j)/(j-cutoff))^2)
    end
    @inbounds @simd for j in 2cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2 * jac
    end
end


const FourierSpectralViscosityCoefficients = Union{Tadmor1989,MadayTadmor1989,TadmorWaagan2012Standard,TadmorWaagan2012Convergent}

function get_parameters(source_of_coefficients::FourierSpectralViscosityCoefficients,
                        D::FourierDerivativeOperator;
                        strength::Real=eltype(D)(1)/size(D,2),
                        cutoff::Int=1+round(Int, sqrt(size(D,2))), #+1: 1 based indexing
                        kwargs...)
    @argcheck cutoff >= 1

    (; strength = strength, cutoff = cutoff)
end


"""
    Tadmor1993()

Coefficients of the Fourier super spectral viscosity given in
- Tadmor (1993)
  Super Viscosity and Spectral Approximations of Nonlinear Conservation Laws.
  Numerical Methods for Fluid Dynamics IV, pp. 69-82.
"""
struct Tadmor1993 <: SourceOfCoefficients end

function Base.show(io::IO, source::Tadmor1993)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
            "Tadmor (1993) \n",
            "  Super Viscosity and Spectral Approximations of Nonlinear Conservation Laws. \n",
            "  Numerical Methods for Fluid Dynamics IV, pp. 69-82.")
    end
end

function set_filter_coefficients!(coefficients::AbstractVector{T},
                                  jac::T, N::Int,
                                  parameters, source::Tadmor1993) where {T<:Real}
    @unpack order, strength, cutoff = parameters
    @argcheck order >= 1
    @argcheck cutoff >= 1
    fill!(coefficients, zero(T))
    jac = jac^2order / N # ^2order: order-th derivative; /N: brfft instead of irfft
    @inbounds @simd for j in cutoff:length(coefficients)
        coefficients[j] = -strength * (j-1)^2order * jac * (1 - (cutoff/j)^2order)
    end
end


const FourierSuperSpectralViscosityCoefficients = Tadmor1993

function get_parameters(source::FourierSuperSpectralViscosityCoefficients,
                        D::FourierDerivativeOperator;
                        order::Int=1,
                        strength=eltype(D)(1)/size(D,2)^(2order-1),
                        cutoff::Int=1+round(Int, size(D,2)^(1-1/2order)), #+1: 1 based indexing
                        kwargs...)
    @argcheck order >= 1
    @argcheck cutoff >= 1

    (; order = order, strength = strength, cutoff = cutoff)
end
