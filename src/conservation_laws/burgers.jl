
"""
    BurgersPeriodicSemidiscretisation

A semidiscretisation of Burgers' equation
    \$\\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0\$
with periodic boundary conditions.
"""
struct BurgersPeriodicSemidiscretisation{T,Derivative<:AbstractDerivativeOperator{T},
                                            Dissipation,
                                            SplitForm<:Union{Val{false}, Val{true}}} <: AbstractSemidiscretisation
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm

    function BurgersPeriodicSemidiscretisation(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm=Val{false}()) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}}
        if dissipation != nothing
            @argcheck size(derivative) == size(dissipation) DimensionMismatch
            @argcheck grid(derivative) == grid(dissipation) ArgumentError
        end
        N = size(derivative, 2)
        tmp1 = Array{T}(undef, N)
        tmp2 = Array{T}(undef, N)
        new{T,Derivative,Dissipation,SplitForm}(derivative, dissipation, tmp1, tmp2, split_form)
    end
end


function Base.show(io::IO, semidisc::BurgersPeriodicSemidiscretisation)
    print(io, "Semidiscretisation of Burgers' equation\n")
    print(io, "  \$ \\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0 \$ \n")
    print(io, "with periodic boundaries using")
    if semidisc.split_form == Val{true}()
        print(io, " a split form and: \n")
    else
        print(io, " no split form and: \n")
    end
    print(io, semidisc.derivative)
    print(io, semidisc.dissipation)
end


function (disc::BurgersPeriodicSemidiscretisation)(du, u, p, t)
    @unpack tmp1, tmp2, derivative, dissipation, split_form = disc
    @boundscheck begin
        @argcheck length(u) == length(tmp1)
        @argcheck length(u) == length(tmp2)
        @argcheck length(u) == length(du)
    end

    # volume terms
    if typeof(split_form) <: Val{true}
        m1_3 = -one(eltype(u)) / 3

        ## u * D * u
        mul!(tmp1, derivative, u)
        @. du = m1_3 * u * tmp1
        ## D * u^2
        @. tmp2 = u^2
        mul!(tmp1, derivative, tmp2)
        @. du += m1_3 * tmp1
    else
        m1_2 = -one(eltype(u)) / 2
        @. tmp2 = u^2
        mul!(tmp1, derivative, tmp2)
        @. du = m1_2 * tmp1
    end

    # dissipation
    if dissipation != nothing
        mul!(tmp1, dissipation, u)
        @. du += tmp1
    end

    nothing
end




"""
    BurgersNonperiodicSemidiscretisation

A semidiscretisation of Burgers' equation
    \$\\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0\$
with boundary conditions `left_bc(t)`, `right_bc(t)`.
"""
struct BurgersNonperiodicSemidiscretisation{T,Derivative<:AbstractDerivativeOperator{T},
                                            Dissipation,
                                            SplitForm<:Union{Val{false}, Val{true}},
                                            LeftBC, RightBC} <: AbstractSemidiscretisation
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm
    left_bc::LeftBC
    right_bc::RightBC

    function BurgersNonperiodicSemidiscretisation(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm, left_bc::LeftBC, right_bc::RightBC) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}, LeftBC, RightBC}
        if dissipation != nothing
            @argcheck size(derivative) == size(dissipation) DimensionMismatch
            @argcheck grid(derivative) == grid(dissipation) ArgumentError
        end
        N = size(derivative, 2)
        tmp1 = Array{T}(undef, N)
        tmp2 = Array{T}(undef, N)
        new{T,Derivative,Dissipation,SplitForm,LeftBC,RightBC}(derivative, dissipation, tmp1, tmp2, split_form, left_bc, right_bc)
    end
end


function Base.show(io::IO, semidisc::BurgersNonperiodicSemidiscretisation)
    print(io, "Semidiscretisation of Burgers' equation\n")
    print(io, "  \$ \\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0 \$ \n")
    print(io, "with nonperiodic boundaries using")
    if semidisc.split_form == Val{true}()
        print(io, " a split form and: \n")
    else
        print(io, " no split form and: \n")
    end
    print(io, semidisc.derivative)
    print(io, semidisc.dissipation)
end


function godunov_flux_burgers(uₗ::T, uᵣ::T) where {T<:Real}
    if uₗ < uᵣ
        if uₗ < 0 && 0 < uᵣ
            zero(T)
        else
            min(uₗ^2/2, uᵣ^2/2)
        end
    else
        max(uₗ^2/2, uᵣ^2/2)
    end
end

function (disc::BurgersNonperiodicSemidiscretisation)(du, u, p, t)
    @unpack tmp1, tmp2, derivative, dissipation, split_form, left_bc, right_bc = disc
    @boundscheck begin
        @argcheck length(u) == length(tmp1)
        @argcheck length(u) == length(tmp2)
        @argcheck length(u) == length(du)
    end

    # volume terms
    if typeof(split_form) <: Val{true}
        m1_3 = -one(eltype(u)) / 3

        ## u * D * u
        mul!(tmp1, derivative, u)
        @. du = m1_3 * u * tmp1
        ## D * u^2
        @. tmp2 = u^2
        mul!(tmp1, derivative, tmp2)
        @. du += m1_3 * tmp1
    else
        m1_2 = -one(eltype(u)) / 2
        @. tmp2 = u^2
        mul!(tmp1, derivative, tmp2)
        @. du = m1_2 * tmp1
    end

    # dissipation
    if dissipation != nothing
        mul!(tmp1, dissipation, u)
        @. du += tmp1
    end

    # boundary conditions via Godunov's flux
    @inbounds fnum_left = godunov_flux_burgers(left_bc(t), u[1])
    @inbounds du[1] += (fnum_left - u[1]^2/2) / left_boundary_weight(derivative)
    @inbounds fnum_right = godunov_flux_burgers(u[end], right_bc(t))
    @inbounds du[end] -= (fnum_right - u[end]^2/2) / right_boundary_weight(derivative)

    nothing
end
