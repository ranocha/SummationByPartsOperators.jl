
"""
    PeriodicDerivativeCoefficients

The coefficients of a derivative operator on a periodic grid.
"""
struct PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients} <: AbstractDerivativeCoefficients{T}
    # coefficients defining the operator and its action
    lower_coef::SVector{LowerOffset, T}
    central_coef::T
    upper_coef::SVector{UpperOffset, T}
    parallel::Parallel
    # corresponding orders etc.
    derivative_order::Int
    accuracy_order  ::Int
    symmetric       ::Bool
    source_of_coefficients::SourceOfCoefficients

    function PeriodicDerivativeCoefficients(lower_coef::SVector{LowerOffset, T}, central_coef::T, upper_coef::SVector{UpperOffset, T},
                                            parallel::Parallel, derivative_order::Int, accuracy_order::Int,
                                            source_of_coefficients::SourceOfCoefficients) where {T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}
        symmetric = LowerOffset == UpperOffset
        if symmetric
            @inbounds for i in Base.OneTo(LowerOffset)
                symmetric = symmetric && lower_coef[i] == upper_coef[i]
            end
        end
        new{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}(lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order, symmetric)
    end
end


# Compute `α*D*u + β*dest` and store the result in `dest`.
function mul!(dest::AbstractVector, coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset}, u::AbstractVector, α, β) where {T,LowerOffset,UpperOffset}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
    end

    @unpack lower_coef, central_coef, upper_coef, parallel = coefficients
    convolve_periodic_boundary_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, β, parallel)
end

# Compute `α*D*u` and store the result in `dest`.
function mul!(dest::AbstractVector, coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset}, u::AbstractVector, α) where {T,LowerOffset,UpperOffset}
    @boundscheck begin
        @argcheck length(dest) == length(u) DimensionMismatch
        @argcheck length(u) > LowerOffset+UpperOffset DimensionMismatch
    end

    @unpack lower_coef, central_coef, upper_coef, parallel = coefficients
    convolve_periodic_boundary_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α)
    convolve_interior_coefficients!(dest, lower_coef, central_coef, upper_coef, u, α, parallel)
end


@generated function convolve_periodic_boundary_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset},
                                                             u::AbstractVector, α, β) where {LowerOffset, UpperOffset}
    ex_lower = :( nothing )
    for i in 1:LowerOffset
        ex = :( lower_coef[$LowerOffset]*u[end-$(LowerOffset-i)] )
        for j in LowerOffset-1:-1:1
            if i-j < 1
                ex = :( $ex + lower_coef[$j]*u[end-$(j-i)] )
            else
                ex = :( $(ex) + lower_coef[$j]*u[$(i-j)] )
            end
        end
        ex = :( $ex + central_coef*u[$i] )
        for j in 1:UpperOffset
            ex = :( $ex + upper_coef[$j]*u[$(i+j)] )
        end
        ex_lower = quote
            $ex_lower
            @inbounds dest[$i] = β*dest[$i] + α*$ex
        end
    end

    ex_upper = :( nothing )
    for i in (UpperOffset-1):-1:0
        ex = :( lower_coef[$LowerOffset]*u[end-$(i+LowerOffset)] )
        for j in LowerOffset-1:-1:1
            ex = :( $ex + lower_coef[$j]*u[end-$(j+i)] )
        end
        ex = :( $ex + central_coef*u[end-$i] )
        for j in 1:UpperOffset
            if i-j < 0
                ex = :( $ex + upper_coef[$j]*u[$(j-i)] )
            else
                ex = :( $ex + upper_coef[$j]*u[end-$(i-j)] )
            end
        end
        ex_upper = quote
            $ex_upper
            @inbounds dest[end-$i] = β*dest[end-$i] + α*$ex
        end
    end

    quote
        $ex_lower
        $ex_upper
    end
end

@generated function convolve_periodic_boundary_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset},
                                                             u::AbstractVector, α) where {LowerOffset, UpperOffset}
    ex_lower = :( nothing )
    for i in 1:LowerOffset
        ex = :( lower_coef[$LowerOffset]*u[end-$(LowerOffset-i)] )
        for j in LowerOffset-1:-1:1
            if i-j < 1
                ex = :( $ex + lower_coef[$j]*u[end-$(j-i)] )
            else
                ex = :( $(ex) + lower_coef[$j]*u[$(i-j)] )
            end
        end
        ex = :( $ex + central_coef*u[$i] )
        for j in 1:UpperOffset
            ex = :( $ex + upper_coef[$j]*u[$(i+j)] )
        end
        ex_lower = quote
            $ex_lower
            @inbounds dest[$i] = α*$ex
        end
    end

    ex_upper = :( nothing )
    for i in (UpperOffset-1):-1:0
        if LowerOffset > 0
            ex = :( lower_coef[$LowerOffset]*u[end-$(i+LowerOffset)] )
            for j in LowerOffset-1:-1:1
                ex = :( $ex + lower_coef[$j]*u[end-$(j+i)] )
            end
            ex = :( $ex + central_coef*u[end-$i] )
        else
            ex = :( central_coef*u[end-$i] )
        end
        for j in 1:UpperOffset
            if i-j < 0
                ex = :( $ex + upper_coef[$j]*u[$(j-i)] )
            else
                ex = :( $ex + upper_coef[$j]*u[end-$(i-j)] )
            end
        end
        ex_upper = quote
            $ex_upper
            @inbounds dest[end-$i] = α*$ex
        end
    end

    quote
        $ex_lower
        $ex_upper
    end
end


@generated function convolve_interior_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset},
                                                    u::AbstractVector, α, β, parallel) where {LowerOffset, UpperOffset}
    ex = :( lower_coef[$LowerOffset]*u[i-$LowerOffset] )
    for j in LowerOffset-1:-1:1
        ex = :( $ex + lower_coef[$j]*u[i-$j] )
    end
    ex = :( $ex + central_coef*u[i] )
    for j in 1:UpperOffset
        ex = :( $ex + upper_coef[$j]*u[i+$j] )
    end

    if parallel <: Val{:threads}
        quote
            Base.@_inline_meta
            @tturbo for i in $(LowerOffset+1):(length(dest)-$UpperOffset)
                dest[i] = β*dest[i] + α*$ex
            end
        end
    else
        quote
            Base.@_inline_meta
            @turbo for i in $(LowerOffset+1):(length(dest)-$UpperOffset)
                dest[i] = β*dest[i] + α*$ex
            end
        end
    end
end

@generated function convolve_interior_coefficients!(dest::AbstractVector, lower_coef::SVector{LowerOffset}, central_coef, upper_coef::SVector{UpperOffset},
                                                    u::AbstractVector, α, parallel) where {LowerOffset, UpperOffset}
    if LowerOffset > 0
        ex = :( lower_coef[$LowerOffset]*u[i-$LowerOffset] )
        for j in LowerOffset-1:-1:1
            ex = :( $ex + lower_coef[$j]*u[i-$j] )
        end
        ex = :( $ex + central_coef*u[i] )
    else
        ex = :( central_coef*u[i] )
    end
    for j in 1:UpperOffset
        ex = :( $ex + upper_coef[$j]*u[i+$j] )
    end

    if parallel <: Val{:threads}
        quote
            Base.@_inline_meta
            @tturbo for i in $(LowerOffset+1):(length(dest)-$UpperOffset)
                dest[i] = α*$ex
            end
        end
    else
        quote
            Base.@_inline_meta
            @turbo for i in $(LowerOffset+1):(length(dest)-$UpperOffset)
                dest[i] = α*$ex
            end
        end
    end
end


"""
    BeljaddLeFlochMishraParés2017()

Coefficients of the periodic operators given in
  Beljadid, LeFloch, Mishra, Parés (2017)
  Schemes with Well-Controlled Dissipation. Hyperbolic Systems in
    Nonconservative Form.
  Communications in Computational Physics 21.4, pp. 913-946.
"""
struct BeljaddLeFlochMishraParés2017 <: SourceOfCoefficients end

function Base.show(io::IO, ::BeljaddLeFlochMishraParés2017)
    print(io,
        "  Beljadid, LeFloch, Mishra, Parés (2017) \n",
        "  Schemes with Well-Controlled Dissipation. Hyperbolic Systems in \n",
        "    Nonconservative Form. \n",
        "  Communications in Computational Physics 21.4, pp. 913-946. \n")
end

"""
    periodic_central_derivative_coefficients(derivative_order, accuracy_order, T=Float64, parallel=Val{:serial}())

Create the `PeriodicDerivativeCoefficients` approximating the `derivative_order`-th
derivative with an order of accuracy `accuracy_order` and scalar type `T`.
The evaluation of the derivative can be parallized using threads by chosing
`parallel=Val{:threads}())`.
"""
function periodic_central_derivative_coefficients(derivative_order, accuracy_order, T=Float64, parallel=Val{:serial}())
    @argcheck derivative_order > 0
    @argcheck accuracy_order > 0
    @argcheck typeof(parallel) <: Union{Val{:serial}, Val{:threads}}

    if isodd(accuracy_order)
        @warn("Central derivative operators support only even orders of accuracy.")
    end

    if derivative_order == 1
        # Exact evaluation of the coefficients, see
        # Beljadid, LeFloch, Mishra, Parés (2017)
        # Schemes with Well-Controlled Dissipation. Hyperbolic Systems in Nonconservative Form.
        # Communications in Computational Physics 21.4, pp. 913-946.
        n = (accuracy_order+1) ÷ 2
        coef = Vector{T}(undef, n)
        for j in 1:n
            tmp = one(Rational{Int128})
            if iseven(j)
                tmp *= -1
            end
            for i in 1:j
                tmp *= (n+1-i) // (n+i)
            end
            coef[j] =  tmp / j
        end
        upper_coef = SVector{n,T}(coef)
        central_coef = zero(T)
        lower_coef = -upper_coef
        source = BeljaddLeFlochMishraParés2017()
    elseif derivative_order == 2
        # Exact evaluation of the coefficients, see
        # Beljadid, LeFloch, Mishra, Parés (2017)
        # Schemes with Well-Controlled Dissipation. Hyperbolic Systems in Nonconservative Form.
        # Communications in Computational Physics 21.4, pp. 913-946.
        n = (accuracy_order+1) ÷ 2
        coef = Vector{T}(undef, n)
        for j in 1:n
            tmp = one(Rational{Int128})
            if iseven(j)
                tmp *= -1
            end
            for i in 1:j
                tmp *= (n+1-i) // (n+i)
            end
            coef[j] =  2*tmp / j^2
        end
        upper_coef = SVector{n,T}(coef)
        central_coef = -2*sum(upper_coef)
        lower_coef = upper_coef
        source = BeljaddLeFlochMishraParés2017()
    elseif derivative_order == 3
        # Exact evaluation of the coefficients, see
        # Beljadid, LeFloch, Mishra, Parés (2017)
        # Schemes with Well-Controlled Dissipation. Hyperbolic Systems in Nonconservative Form.
        # Communications in Computational Physics 21.4, pp. 913-946.
        n = (accuracy_order+1) ÷ 2
        coef = Vector{T}(undef, n)
        for j in 1:n
            tmp = one(Rational{Int128})
            if isodd(j)
                tmp *= -1
            end
            for i in 1:j
                tmp *= (n+1-i) // (n+i)
            end
            fac = zero(Rational{Int128})
            for k in 1:n
                k == j && continue
                fac += 1 // k^2
            end
            coef[j] =  6*tmp*fac/ j
        end
        upper_coef = SVector{n,T}(coef)
        central_coef = zero(T)
        lower_coef = -upper_coef
        source = BeljaddLeFlochMishraParés2017()
    else
        throw(ArgumentError("Derivative order $derivative_order not implemented yet."))
    end

    PeriodicDerivativeCoefficients(lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order, source)
end

"""
    Fornberg1998()

Coefficients of the periodic operators given in
  Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.
"""
struct Fornberg1998 <: SourceOfCoefficients end

function Base.show(io::IO, ::Fornberg1998)
    print(io,
        "  Fornberg (1998) \n",
        "  Calculation of Weights in Finite Difference Formulas. \n",
        "  SIAM Rev. 40.3, pp. 685-691. \n")
end

"""
    fornberg(x::Vector{T}, m::Int) where {T}

Calculate the weights of a finite difference approximation of the `m`th derivative
with maximal order of accuracy at `0` using the nodes `x`, see
Fornberg (1998)
Calculation of Weights in Finite Difference Formulas
SIAM Rev. 40.3, pp. 685-691.
"""
function fornberg(xx::Vector{T}, m::Int) where {T}
    x = sort(xx)
    @argcheck x[1] <= 0
    @argcheck x[end] >= 0
    @argcheck length(unique(x)) == length(x)
    @argcheck zero(T) in x
    @argcheck m >= 1

    z = zero(T)
    n = length(x) - 1
    c = fill(zero(T), length(x), m+1)
    c1 = one(T)
    c4 = x[1] - z
    c[1,1] = one(T)
    for i in 1:n
        mn = min(i,m)
        c2 = one(T)
        c5 = c4
        c4 = x[i+1] - z
        for j in 0:i-1
            c3 = x[i+1] - x[j+1]
            c2 = c2 * c3
            if j == i-1
                for k in mn:-1:1
                    c[i+1,k+1] = c1 * (k*c[i,k]-c5*c[i,k+1]) / c2
                end
                c[i+1,1] = -c1*c5*c[i,1] / c2
            end
            for k in mn:-1:1
                c[j+1,k+1] = (c4*c[j+1,k+1]-k*c[j+1,k]) / c3
            end
            c[j+1,1] = c4*c[j+1,1] / c3
        end
        c1 = c2
    end

    c[:,end]
end

"""
    periodic_derivative_coefficients(derivative_order, accuracy_order, left_offset=-(accuracy_order+1)÷2, T=Float64, parallel=Val{:serial}())

Create the `PeriodicDerivativeCoefficients` approximating the `derivative_order`-th
derivative with an order of accuracy `accuracy_order` and scalar type `T` where
the leftmost grid point used is determined by `left_offset`.
The evaluation of the derivative can be parallized using threads by chosing
parallel=Val{:threads}())`.
"""
function periodic_derivative_coefficients(derivative_order, accuracy_order, left_offset::Int=-(accuracy_order+1)÷2, T=Float64, parallel=Val{:serial}())
    @argcheck -accuracy_order <= left_offset <= 0

    x = Rational{Int128}.(left_offset:left_offset+accuracy_order)
    c = fornberg(x, derivative_order)

    z = zero(eltype(x))
    lower_idx = x .< z
    central_idx = first(findall(xx->xx==z, x))
    upper_idx = x .> z

    LowerOffset = sum(lower_idx)
    UpperOffset = sum(upper_idx)

    lower_coef = SVector{LowerOffset, T}(T.(reverse(c[lower_idx])))
    central_coef = T(c[central_idx])
    upper_coef = SVector{UpperOffset, T}(T.(c[upper_idx]))

    source = Fornberg1998()

    PeriodicDerivativeCoefficients(lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order, source)
end



"""
    Holoborodko2008()

Coefficients of the periodic operators given in
  Holoborodko (2008)
  Smooth Noise Robust Differentiators.
  http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
"""
struct Holoborodko2008 <: SourceOfCoefficients end

function Base.show(io::IO, ::Holoborodko2008)
    print(io,
        "  Holoborodko (2008) \n",
        "  Smooth Noise Robust Differentiators. \n",
        "  http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/ \n")
end

"""
    periodic_derivative_coefficients(source::Holoborodko2008, derivative_order, accuracy_order;
                                     T=Float64, parallel=Val{:serial}(),
                                     stencil_width=accuracy_order+3)

Create the `PeriodicDerivativeCoefficients` approximating the `derivative_order`-th
derivative with an order of accuracy `accuracy_order` and scalar type `T` given
by Holoborodko2008.
The evaluation of the derivative can be parallized using threads by chosing
parallel=Val{:threads}())`.
"""
function periodic_derivative_coefficients(source::Holoborodko2008, derivative_order, accuracy_order;
                                          T=Float64, parallel=Val{:serial}(),
                                          stencil_width=accuracy_order+3)
    method_exists = true
    if derivative_order == 1
        if accuracy_order == 2
            if stencil_width == 5
                lower_coef = SVector{2, T}([
                    -2//8, -1//8
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            elseif stencil_width == 7
                lower_coef = SVector{3, T}([
                    -5//32, -4//32, -1//32
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            elseif stencil_width == 9
                lower_coef = SVector{4, T}([
                    -14//128, -14//128, -6//128, -1//128
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            elseif stencil_width == 11
                lower_coef = SVector{5, T}([
                    -42//512, -48//512, -27//512, -8//512, -1//512
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            else
                method_exists = false
            end
        elseif accuracy_order == 4
            if stencil_width == 7
                lower_coef = SVector{3, T}([
                    -39//96, -12//96, 5//96
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            elseif stencil_width == 9
                lower_coef = SVector{4, T}([
                    -27//96, -16//96, 1//96, 2//96
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            elseif stencil_width == 11
                lower_coef = SVector{5, T}([
                    -322//1536, -256//1536, -39//1536, 32//1536, 11//1536
                ])
                upper_coef = -lower_coef
                central_coef = T(0)
            else
                method_exists = false
            end
        else
            method_exists = false
        end
    elseif derivative_order == 2
        if accuracy_order == 2
            if stencil_width == 5
                lower_coef = SVector{2, T}([
                    0, 1//4
                ])
                upper_coef = lower_coef
                central_coef = T(-2//4)
            elseif stencil_width == 7
                lower_coef = SVector{3, T}([
                    -1//16, 2//16, 1//16
                ])
                upper_coef = lower_coef
                central_coef = T(-4//16)
            elseif stencil_width == 9
                lower_coef = SVector{4, T}([
                    -4//64, 4//64, 4//64, 1//64
                ])
                upper_coef = lower_coef
                central_coef = T(-10//64)
            else
                method_exists = false
            end
        elseif accuracy_order == 4
            if stencil_width == 7
                lower_coef = SVector{3, T}([
                    1//12, 5//12, -1//12
                ])
                upper_coef = lower_coef
                central_coef = T(-10//12)
            elseif stencil_width == 9
                lower_coef = SVector{4, T}([
                    -12//192, 52//192, 12//192, -7//192
                ])
                upper_coef = lower_coef
                central_coef = T(-90//192)
            else
                method_exists = false
            end
        else
            method_exists = false
        end
    else
        method_exists = false
    end
    if method_exists == false
        throw(ArgumentError("Method with derivative_order=$derivative_order, accuracy_order=$accuracy_order, stencil_width=$stencil_width not implemented/derived."))
    end

    PeriodicDerivativeCoefficients(lower_coef, central_coef, upper_coef, parallel, derivative_order, accuracy_order, source)
end


"""
    PeriodicDerivativeOperator

A derivative operator on a uniform periodic grid.
"""
struct PeriodicDerivativeOperator{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid} <: AbstractPeriodicDerivativeOperator{T}
    coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients}
    grid_compute::Grid   # N   nodes, including the left and excluding the right boundary
    grid_evaluate::Grid #  N+1 nodes, including both boundaries
    Δx::T
    factor::T

    function PeriodicDerivativeOperator(coefficients::PeriodicDerivativeCoefficients{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients},
                                        grid_evaluate::Grid) where {T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}
        @argcheck length(grid_evaluate) > LowerOffset+UpperOffset DimensionMismatch
        grid_compute = range(grid_evaluate[1], stop=grid_evaluate[end-1], length=length(grid_evaluate)-1)

        Δx = step(grid_evaluate)
        factor = inv(Δx^coefficients.derivative_order)
        new{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}(coefficients, grid_compute, grid_evaluate, Δx, factor)
    end
end

@inline source_of_coefficients(D::PeriodicDerivativeOperator) = source_of_coefficients(D.coefficients)


function Base.show(io::IO, D::PeriodicDerivativeOperator{T,LowerOffset,UpperOffset}) where {T,LowerOffset,UpperOffset}
    x = D.grid_evaluate
    if derivative_order(D) == 1
        print(io, "Periodic 1st derivative operator of order ")
    elseif  derivative_order(D) == 2
        print(io, "Periodic 2nd derivative operator of order ")
    elseif  derivative_order(D) == 3
        print(io, "Periodic 3rd derivative operator of order ")
    else
        print(io, "Periodic ", derivative_order(D), "th derivative operator of order ")
    end
    print(io, accuracy_order(D), " {T=", T, ", Parallel=", typeof(D.coefficients.parallel), "} \n")
    print(io, "on a grid in [", first(x), ", ", last(x), "] using ", size(D, 1), " nodes, \n")
    print(io, "stencils with ", LowerOffset, " nodes to the left, ", UpperOffset,
                " nodes to the right, and coefficients from \n", source_of_coefficients(D))
end


# Compute `α*D*u + β*dest` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::PeriodicDerivativeOperator, u::AbstractVector, α, β)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, α*D.factor, β)
end

# Compute `α*D*u` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, D::PeriodicDerivativeOperator, u::AbstractVector, α)
    @boundscheck begin
        @argcheck size(D, 2) == length(u) DimensionMismatch
        @argcheck size(D, 1) == length(dest) DimensionMismatch
    end
    @inbounds mul!(dest, D.coefficients, u, α*D.factor)
end


"""
    integrate(func, u, D::PeriodicDerivativeOperator)

Map the function `func` to the coefficients `u` and integrate with respect to
the quadrature rule associated with the periodic derivative operator `D`.
"""
function integrate(func, u::AbstractVector, D::PeriodicDerivativeOperator)
    @boundscheck begin
        length(u) == length(grid(D))
    end
    @unpack Δx = D

    @inbounds res = sum(func, u)

    Δx * res
end

function mass_matrix(D::PeriodicDerivativeOperator)
    @unpack Δx = D

    Δx * I
end


"""
    periodic_central_derivative_operator(derivative_order, accuracy_order,
                                         xmin, xmax, N, parallel=Val{:serial}())

Create a [`PeriodicDerivativeOperator`](@ref) approximating the `derivative_order`-th
derivative on a uniform grid between `xmin` and `xmax` with `N` grid points up
to order of accuracy `accuracy_order`.
The evaluation of the derivative can be parallized using threads by chosing
`parallel=Val{:threads}())`.
"""
function periodic_central_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N, parallel=Val{:serial}())
    grid_evaluate = range(xmin, stop=xmax, length=N+1) # N includes two identical boundary nodes
    coefficients = periodic_central_derivative_coefficients(derivative_order, accuracy_order, eltype(grid_evaluate), parallel)
    PeriodicDerivativeOperator(coefficients, grid_evaluate)
end

"""
    periodic_derivative_operator(derivative_order, accuracy_order,
                                 xmin, xmax, N,
                                 left_offset=-(accuracy_order+1)÷2,
                                 parallel=Val{:serial}())
    periodic_derivative_operator(; derivative_order, accuracy_order,
                                 xmin, xmax, N,
                                 left_offset=-(accuracy_order+1)÷2,
                                 parallel=Val{:serial}())

Create a [`PeriodicDerivativeOperator`](@ref) approximating the `derivative_order`-th
derivative on a uniform grid between `xmin` and `xmax` with `N` grid points up
to order of accuracy `accuracy_order` where the leftmost grid point used is
determined by `left_offset`.
The evaluation of the derivative can be parallized using threads by chosing
`parallel=Val{:threads}())`.

## Examples

```jldoctest
julia> periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                    xmin=0.0, xmax=1.0, N=11)
Periodic 1st derivative operator of order 2 {T=Float64, Parallel=Val{:serial}}
on a grid in [0.0, 1.0] using 11 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients from
  Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

```
"""
function periodic_derivative_operator(derivative_order::Integer, accuracy_order,
                                      xmin, xmax, N, left_offset::Int=-(accuracy_order+1)÷2,
                                      parallel::Union{Val{:serial},Val{:threads}}=Val{:serial}())
    grid = range(xmin, stop=xmax, length=N+1) # N includes two identical boundary nodes
    coefficients = periodic_derivative_coefficients(derivative_order, accuracy_order, left_offset, eltype(grid), parallel)
    PeriodicDerivativeOperator(coefficients, grid)
end

@inline function periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N,
                                              parallel::Union{Val{:serial},Val{:threads}}, left_offset::Int=-(accuracy_order+1)÷2)
    periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N, left_offset, parallel)
end

function periodic_derivative_operator(; derivative_order, accuracy_order,
                                      xmin, xmax, N, left_offset::Int=-(accuracy_order+1)÷2,
                                      parallel::Union{Val{:serial},Val{:threads}}=Val{:serial}())
    periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N,
                                 left_offset, parallel)
end

"""
    periodic_derivative_operator(source::Holoborodko2008,
                                 derivative_order, accuracy_order,
                                 xmin, xmax, N; parallel=Val{:serial}(), kwargs...)
    periodic_derivative_operator(source::Holoborodko2008;
                                 derivative_order, accuracy_order,
                                 xmin, xmax, N, parallel=Val{:serial}(), kwargs...)

Create a `PeriodicDerivativeOperator` approximating the `derivative_order`-th
derivative on a uniform grid between `xmin` and `xmax` with `N` grid points up
to order of accuracy `accuracy_order` where the leftmost grid point used is
determined by `left_offset`.
The evaluation of the derivative can be parallized using threads by chosing
`parallel=Val{:threads}())`.

## Examples

```jldoctest
julia> periodic_derivative_operator(Holoborodko2008(), derivative_order=1, accuracy_order=2,
                                    xmin=0.0, xmax=1.0, N=11)
Periodic 1st derivative operator of order 2 {T=Float64, Parallel=Val{:serial}}
on a grid in [0.0, 1.0] using 11 nodes,
stencils with 2 nodes to the left, 2 nodes to the right, and coefficients from
  Holoborodko (2008)
  Smooth Noise Robust Differentiators.
  http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/

```
"""
function periodic_derivative_operator(source::Holoborodko2008, derivative_order, accuracy_order,
                                      xmin, xmax, N; parallel=Val{:serial}(), kwargs...)
    grid = range(xmin, stop=xmax, length=N+1) # N includes two identical boundary nodes
    coefficients = periodic_derivative_coefficients(source, derivative_order, accuracy_order;
                                                    kwargs..., T=eltype(grid), parallel=parallel)
    PeriodicDerivativeOperator(coefficients, grid)
end

function periodic_derivative_operator(source::Holoborodko2008; derivative_order, accuracy_order,
                                      xmin, xmax, N, parallel=Val{:serial}(), kwargs...)
    periodic_derivative_operator(source, derivative_order, accuracy_order, xmin, xmax, N;
                                 parallel=parallel, kwargs...)
end

"""
    periodic_derivative_operator(derivative_order, accuracy_order, grid,
                                 left_offset=-(accuracy_order+1)÷2, parallel=Val{:serial}())

Create a `PeriodicDerivativeOperator` approximating the `derivative_order`-th
derivative on thr uniform `grid` up to order of accuracy `accuracy_order` where
the leftmost grid point used is determined by `left_offset`.
The evaluation of the derivative can be parallized using threads by chosing
`parallel=Val{:threads}())`.
"""
function periodic_derivative_operator(derivative_order, accuracy_order, grid::Union{LinRange,StepRangeLen},
                                      left_offset::Int=-(accuracy_order+1)÷2, parallel=Val{:serial}())
    coefficients = periodic_derivative_coefficients(derivative_order, accuracy_order, left_offset, eltype(grid), parallel)
    PeriodicDerivativeOperator(coefficients, grid)
end

@inline function periodic_derivative_operator(derivative_order, accuracy_order, grid::Union{LinRange,StepRangeLen},
                                              parallel::Union{Val{:serial},Val{:threads}}, left_offset::Int=-(accuracy_order+1)÷2)
    periodic_derivative_operator(derivative_order, accuracy_order, grid, left_offset, parallel)
end



"""
    PeriodicDissipationOperator

A dissipation operator on a periodic finite difference grid.
"""
struct PeriodicDissipationOperator{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid} <: AbstractPeriodicDerivativeOperator{T}
    factor::T
    Di::PeriodicDerivativeOperator{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}

    function PeriodicDissipationOperator(factor::T, Di::PeriodicDerivativeOperator{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}) where {T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}
        # check validity
        @argcheck iseven(derivative_order(Di))
        @argcheck (-1)^(1+derivative_order(Di)÷2)*factor >= 0
        new{T,LowerOffset,UpperOffset,Parallel,SourceOfCoefficients,Grid}(factor, Di)
    end
end

grid(Di::PeriodicDissipationOperator) = grid(Di.Di)
derivative_order(Di::PeriodicDissipationOperator) = derivative_order(Di.Di)
accuracy_order(Di::PeriodicDissipationOperator) = accuracy_order(Di.Di)
source_of_coefficients(Di::PeriodicDissipationOperator) = MattssonSvärdNordström2004()

function Base.show(io::IO, Di::PeriodicDissipationOperator{T}) where {T}
    if  derivative_order(Di) == 2
        print(io, "SBP 2nd derivative dissipation operator of order ")
    else
        print(io, "SBP ", derivative_order(Di), "th derivative dissipation operator of order ")
    end
    x = Di.Di.grid_evaluate
    print(io, accuracy_order(Di), " {T=", T, ", Parallel=", typeof(Di.Di.coefficients.parallel), "} \n")
    print(io, "on a grid in [", first(x), ", ", last(x), "] using ", size(Di, 1), " nodes \n")
    print(io, "and coefficients given in \n")
    print(io, source_of_coefficients(Di))
end


# Compute `α*D*u + β*dest` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, Di::PeriodicDissipationOperator,
                                       u::AbstractVector, α, β)
    mul!(dest, Di.Di, u, Di.factor*α, β)
end

# Compute `α*D*u` and store the result in `dest`.
Base.@propagate_inbounds function mul!(dest::AbstractVector, Di::PeriodicDissipationOperator,
                                       u::AbstractVector, α)
    mul!(dest, Di.Di, u, Di.factor*α)
end


"""
    dissipation_operator(D::PeriodicDerivativeOperator;
                         strength=one(eltype(D)),
                         order=accuracy_order(D),
                         parallel=D.coefficients.parallel)

Create a negative semidefinite `DissipationOperator` using undivided differences
approximating a `order`-th derivative with strength `strength` adapted to the
derivative operator `D`.
The evaluation of the derivative can be parallized using threads by chosing
`parallel=Val{:threads}())`.
"""
function dissipation_operator(D::PeriodicDerivativeOperator;
                              strength=eltype(D)(1),
                              order=accuracy_order(D),
                              parallel=D.coefficients.parallel)
    @argcheck iseven(order) ArgumentError("Dissipation operators require even derivatives.")
    # account for the grid spacing
    factor = strength * (-1)^(1 + order÷2) * D.Δx^order
    x = D.grid_evaluate
    Di = periodic_derivative_operator(order, order, first(x), last(x), size(D, 1), parallel)
    PeriodicDissipationOperator(factor, Di)
end




struct PeriodicRationalDerivativeOperator{T<:Real, Dtype<:PeriodicDerivativeOperator{T}, Nnum, Nden, RFFT, IRFFT} <: AbstractPeriodicDerivativeOperator{T}
    D::Dtype
    num_coef::NTuple{Nnum,T}
    den_coef::NTuple{Nden,T}
    tmp::Vector{Complex{T}}
    eigval::Vector{Complex{T}}
    rfft_plan::RFFT
    irfft_plan::IRFFT

    function PeriodicRationalDerivativeOperator(D::Dtype, num_coef::NTuple{Nnum,T}, den_coef::NTuple{Nden,T}, tmp::Vector{Complex{T}}, eigval::Vector{Complex{T}}, rfft_plan::RFFT, irfft_plan::IRFFT) where {T<:Real, Dtype<:PeriodicDerivativeOperator{T}, Nnum, Nden, RFFT, IRFFT}
        @argcheck length(irfft_plan) == length(tmp) DimensionMismatch
        @argcheck length(irfft_plan) == length(eigval) DimensionMismatch
        @argcheck length(irfft_plan) == (length(rfft_plan)÷2)+1 DimensionMismatch
        @argcheck length(grid(D)) == length(rfft_plan) DimensionMismatch

        new{T,Dtype,Nnum,Nden,RFFT,IRFFT}(D, num_coef, den_coef, tmp, eigval, rfft_plan, irfft_plan)
    end
end


function _eigvals!(eigval::Vector{Complex{T}}, D::PeriodicDerivativeOperator{T}) where {T<:Real}
    @unpack factor = D
    @unpack lower_coef, central_coef, upper_coef = D.coefficients
    N = length(grid(D))

    @inbounds for idx in 1:length(eigval)
        tmp_eigval = zero(T)
        for j in Base.OneTo(length(lower_coef))
            tmp_eigval += factor * lower_coef[j] * (exp(2π * im * (idx-1) * (-j) / N))
        end
        tmp_eigval += factor * central_coef
        for j in Base.OneTo(length(upper_coef))
            tmp_eigval += factor * upper_coef[j] * (exp(2π * im * (idx-1) * (+j) / N))
        end
        eigval[idx] = tmp_eigval
    end

    eigval
end

function PeriodicRationalDerivativeOperator(D::PeriodicDerivativeOperator{T}) where {T<:Real}
    x = grid(D)
    u = zero.(x)
    rfft_plan = plan_rfft(u)
    tmp = rfft_plan * u
    irfft_plan = plan_irfft(tmp, length(x))

    eigval = similar(tmp)
    _eigvals!(eigval, D)

    PeriodicRationalDerivativeOperator(D, (zero(T), one(T)), (one(T),), tmp, eigval, rfft_plan, irfft_plan)
end

function PeriodicRationalDerivativeOperator(rat::PeriodicRationalDerivativeOperator, num_coef, den_coef)

    PeriodicRationalDerivativeOperator(rat.D, num_coef, den_coef, rat.tmp, rat.eigval, rat.rfft_plan, rat.irfft_plan)
end

function PeriodicRationalDerivativeOperator(D::PeriodicDerivativeOperator, num_coef, den_coef)

    PeriodicRationalDerivativeOperator(PeriodicRationalDerivativeOperator(D), num_coef, den_coef)
end

Base.size(rat::PeriodicRationalDerivativeOperator) = size(rat.D)
grid(rat::PeriodicRationalDerivativeOperator) = grid(rat.D)

function Base.show(io::IO, rat::PeriodicRationalDerivativeOperator)
    print(io, "Rational periodic operator with coefficients\n")
    print(io, rat.num_coef)
    print(io, "\nand\n")
    print(io, rat.den_coef)
    print(io, "\nof the operator:\n")
    print(io, rat.D)
end


function mul_poly(coef1, coef2)
    T = promote_type(eltype(coef1), eltype(coef2))

    coef = ntuple(idx->zero(T), (length(coef1)-1) + (length(coef2)-1) + 1)
    for idx1 in 0:length(coef1)-1, idx2 in 0:length(coef2)-1
        coef = Base.setindex(coef, coef[idx1+idx2+1] + coef1[idx1+1] * coef2[idx2+1], idx1+idx2+1)
    end

    coef
end

function add_poly(coef1, coef2)
    T = promote_type(eltype(coef1), eltype(coef2))

    coef = ntuple(idx->zero(T), max(length(coef1), length(coef2)))
    for idx in 1:length(coef1)
        coef = Base.setindex(coef, coef[idx] + coef1[idx], idx)
    end
    for idx in 1:length(coef2)
        coef = Base.setindex(coef, coef[idx] + coef2[idx], idx)
    end

    coef
end

function subtract_poly(coef1, coef2)
    T = promote_type(eltype(coef1), eltype(coef2))

    coef = ntuple(idx->zero(T), max(length(coef1), length(coef2)))
    for idx in 1:length(coef1)
        coef = Base.setindex(coef, coef[idx] + coef1[idx], idx)
    end
    for idx in 1:length(coef2)
        coef = Base.setindex(coef, coef[idx] - coef2[idx], idx)
    end

    coef
end


function Base.inv(rat::PeriodicRationalDerivativeOperator)
    PeriodicRationalDerivativeOperator(rat, rat.den_coef, rat.num_coef)
end

function Base.:*(D1::PeriodicDerivativeOperator, D2::PeriodicDerivativeOperator)
    T = eltype(D1)
    @argcheck D1 == D2 ArgumentError

    PeriodicRationalDerivativeOperator(D1, (zero(T), zero(T), one(T)), (one(T),))
end

function Base.literal_pow(::typeof(^), D1::PeriodicDerivativeOperator, ::Val{P}) where {P}
    T = eltype(D1)
    coef = Base.setindex( ntuple(_->zero(T), Val{P+1}()), one(T), P+1)
    PeriodicRationalDerivativeOperator(D1, coef, (one(T),))
end

function Base.:*(factor::Union{Real,Integer}, rat::PeriodicRationalDerivativeOperator)
    @unpack num_coef = rat
    for idx in 1:length(num_coef)
        num_coef = Base.setindex(num_coef, factor*num_coef[idx], idx)
    end

    PeriodicRationalDerivativeOperator(rat, num_coef, rat.den_coef)
end

function Base.:*(rat::PeriodicRationalDerivativeOperator, factor::Union{Real,Integer})
    factor * rat
end

function Base.:*(D::PeriodicDerivativeOperator, factor::Union{Real,Integer})
    PeriodicRationalDerivativeOperator(D) * factor
end

function Base.:*(factor::Union{Real,Integer}, D::PeriodicDerivativeOperator)
    D * factor
end

function Base.:*(rat::PeriodicRationalDerivativeOperator, scaling::UniformScaling)
    scaling.λ * rat
end

function Base.:*(scaling::UniformScaling, rat::PeriodicRationalDerivativeOperator)
    rat * scaling
end

function Base.:*(D::PeriodicDerivativeOperator, scaling::UniformScaling)
    scaling * PeriodicRationalDerivativeOperator(D)
end

function Base.:*(scaling::UniformScaling, D::PeriodicDerivativeOperator)
    PeriodicRationalDerivativeOperator(D) * scaling
end


function Base.:+(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    @argcheck rat1.D == rat2.D ArgumentError

    num_coef = add_poly(mul_poly(rat1.num_coef, rat2.den_coef), mul_poly(rat1.den_coef, rat2.num_coef))
    den_coef = mul_poly(rat1.den_coef, rat2.den_coef)
    PeriodicRationalDerivativeOperator(rat1, num_coef, den_coef)
end

function Base.:+(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicDerivativeOperator)
    rat1 + PeriodicRationalDerivativeOperator(rat2)
end

function Base.:+(rat1::PeriodicDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    PeriodicRationalDerivativeOperator(rat1) + rat2
end

function Base.:+(rat::PeriodicRationalDerivativeOperator, scaling::UniformScaling)
    @unpack num_coef, den_coef = rat
    num_coef = add_poly(num_coef, den_coef)

    PeriodicRationalDerivativeOperator(rat, num_coef, den_coef)
end

function Base.:+(scaling::UniformScaling, rat::PeriodicRationalDerivativeOperator)
    rat + scaling
end

function Base.:+(D::PeriodicDerivativeOperator, scaling::UniformScaling)
    PeriodicRationalDerivativeOperator(D) + scaling
end

function Base.:+(scaling::UniformScaling, D::PeriodicDerivativeOperator)
    scaling + PeriodicRationalDerivativeOperator(D)
end

function Base.:-(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    @argcheck rat1.D == rat2.D ArgumentError

    num_coef = subtract_poly(mul_poly(rat1.num_coef, rat2.den_coef), mul_poly(rat1.den_coef, rat2.num_coef))
    den_coef = mul_poly(rat1.den_coef, rat2.den_coef)
    PeriodicRationalDerivativeOperator(rat1, num_coef, den_coef)
end

function Base.:-(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicDerivativeOperator)
    rat1 - PeriodicRationalDerivativeOperator(rat2)
end

function Base.:-(rat1::PeriodicDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    PeriodicRationalDerivativeOperator(rat1) - rat2
end

function Base.:-(rat::PeriodicRationalDerivativeOperator, scaling::UniformScaling)
    @unpack num_coef, den_coef = rat
    num_coef = subtract_poly(num_coef, den_coef)

    PeriodicRationalDerivativeOperator(rat, num_coef, den_coef)
end

function Base.:-(scaling::UniformScaling, rat::PeriodicRationalDerivativeOperator)
    @unpack num_coef, den_coef = rat
    num_coef = subtract_poly(den_coef, num_coef)

    PeriodicRationalDerivativeOperator(rat, num_coef, den_coef)
end

function Base.:-(D::PeriodicDerivativeOperator, scaling::UniformScaling)
    PeriodicRationalDerivativeOperator(D) - scaling
end

function Base.:-(scaling::UniformScaling, D::PeriodicDerivativeOperator)
    scaling - PeriodicRationalDerivativeOperator(D)
end

function Base.:*(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    @argcheck rat1.D == rat2.D ArgumentError

    num_coef = mul_poly(rat1.num_coef, rat2.num_coef)
    den_coef = mul_poly(rat1.den_coef, rat2.den_coef)
    PeriodicRationalDerivativeOperator(rat1, num_coef, den_coef)
end

function Base.:*(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicDerivativeOperator)
    rat1 * PeriodicRationalDerivativeOperator(rat2)
end

function Base.:*(rat1::PeriodicDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    PeriodicRationalDerivativeOperator(rat1) * rat2
end

function Base.:/(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    @argcheck rat1.D == rat2.D ArgumentError

    num_coef = mul_poly(rat1.num_coef, rat2.den_coef)
    den_coef = mul_poly(rat1.den_coef, rat2.num_coef)
    PeriodicRationalDerivativeOperator(rat1, num_coef, den_coef)
end

function Base.:/(rat1::PeriodicRationalDerivativeOperator, rat2::PeriodicDerivativeOperator)
    rat1 / PeriodicRationalDerivativeOperator(rat2)
end

function Base.:/(rat1::PeriodicDerivativeOperator, rat2::PeriodicRationalDerivativeOperator)
    PeriodicRationalDerivativeOperator(rat1) / rat2
end


function mul!(dest::AbstractVector{T}, rat::PeriodicRationalDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack D, num_coef, den_coef, tmp, eigval, rfft_plan, irfft_plan = rat
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp))
        tmp[j] *= evalpoly(eigval[j], num_coef) / evalpoly(eigval[j], den_coef)
    end
    mul!(dest, irfft_plan, tmp)
end

function LinearAlgebra.ldiv!(dest::AbstractVector{T}, rat::PeriodicRationalDerivativeOperator, u::AbstractVector{T}) where {T}
    @unpack D, num_coef, den_coef, tmp, eigval, rfft_plan, irfft_plan = rat
    N, _ = size(D)
    @boundscheck begin
        @argcheck N == length(u)
        @argcheck N == length(dest)
    end

    mul!(tmp, rfft_plan, u)
    @inbounds @simd for j in Base.OneTo(length(tmp))
        tmp[j] *= evalpoly(eigval[j], den_coef) / evalpoly(eigval[j], num_coef)
    end
    mul!(dest, irfft_plan, tmp)
end

function Base.:\(rat::PeriodicRationalDerivativeOperator{T}, u::AbstractVector{T}) where {T}
    dest = similar(u)
    ldiv!(dest, rat, u)
end
