
"""
    BurgersPeriodicSemidiscretization(D, Di, split_form)

A semidiscretization of Burgers' equation
    \$\\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0\$
with periodic boundary conditions.

`D` is a first-derivative SBP operator, `Di` an associated dissipation operator
or `nothing`, and `split_form::Union{Val(true), Val(false)}` determines whether
the canonical split form or the conservative form is used.
"""
@auto_hash_equals struct BurgersPeriodicSemidiscretization{T, Derivative<:AbstractDerivativeOperator{T},
                                                              Dissipation,
                                                              SplitForm<:Union{Val{false}, Val{true}}} <: AbstractSemidiscretization
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm

    function BurgersPeriodicSemidiscretization(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm=Val{false}()) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}}
        if dissipation !== nothing
            @argcheck size(derivative) == size(dissipation) DimensionMismatch
            @argcheck grid(derivative) == grid(dissipation) ArgumentError
        end
        N = size(derivative, 2)
        tmp1 = Array{T}(undef, N)
        tmp2 = Array{T}(undef, N)
        new{T,Derivative,Dissipation,SplitForm}(derivative, dissipation, tmp1, tmp2, split_form)
    end
end


function Base.show(io::IO, semi::BurgersPeriodicSemidiscretization)
    if get(io, :compact, false)
        print(io, "Semidiscretization of Burgers' equation (periodic)")
    else
        print(io, "Semidiscretization of Burgers' equation\n")
        print(io, "  \$ \\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0 \$ \n")
        print(io, "with periodic boundaries using")
        if semi.split_form == Val{true}()
            print(io, " a split form and: \n")
        else
            print(io, " no split form and: \n")
        end
        print(io, semi.derivative)
        print(io, semi.dissipation)
    end
end


function (disc::BurgersPeriodicSemidiscretization)(du, u, p, t)
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
    if dissipation !== nothing
        mul!(tmp1, dissipation, u)
        @. du += tmp1
    end

    nothing
end




"""
    BurgersNonperiodicSemidiscretization(D, Di, split_form, left_bc, right_bc)

A semidiscretization of Burgers' equation
    \$\\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0\$
with boundary conditions `left_bc(t)`, `right_bc(t)`.

`D` is a first-derivative SBP operator, `Di` an associated dissipation operator
or `nothing`, and `split_form::Union{Val(true), Val(false)}` determines whether
the canonical split form or the conservative form is used.
"""
@auto_hash_equals struct BurgersNonperiodicSemidiscretization{T, Derivative<:AbstractDerivativeOperator{T},
                                                              Dissipation,
                                                              SplitForm<:Union{Val{false}, Val{true}},
                                                              LeftBC, RightBC} <: AbstractSemidiscretization
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm
    left_bc::LeftBC
    right_bc::RightBC

    function BurgersNonperiodicSemidiscretization(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm, left_bc::LeftBC, right_bc::RightBC) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}, LeftBC, RightBC}
        if dissipation !== nothing
            @argcheck size(derivative) == size(dissipation) DimensionMismatch
            @argcheck grid(derivative) == grid(dissipation) ArgumentError
        end
        N = size(derivative, 2)
        tmp1 = Array{T}(undef, N)
        tmp2 = Array{T}(undef, N)
        new{T,Derivative,Dissipation,SplitForm,LeftBC,RightBC}(derivative, dissipation, tmp1, tmp2, split_form, left_bc, right_bc)
    end
end


function Base.show(io::IO, semi::BurgersNonperiodicSemidiscretization)
    if get(io, :compact, false)
        print(io, "Semidiscretization of Burgers' equation (non-periodic)")
    else
        print(io, "Semidiscretization of Burgers' equation\n")
        print(io, "  \$ \\partial_t u(t,x) + \\partial_x \\frac{u(t,x)^2}{2} = 0 \$ \n")
        print(io, "with nonperiodic boundaries using")
        if semi.split_form == Val{true}()
            print(io, " a split form and: \n")
        else
            print(io, " no split form and: \n")
        end
        print(io, semi.derivative)
        print(io, semi.dissipation)
    end
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

function (disc::BurgersNonperiodicSemidiscretization)(du, u, p, t)
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
    if dissipation !== nothing
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
