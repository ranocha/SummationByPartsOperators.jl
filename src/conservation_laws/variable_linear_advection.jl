"""
    VariableLinearAdvectionNonperiodicSemidiscretization(D, Di, a, split_form,
                                                         left_bc, right_bc)

A semidiscretization of the linear advection equation
    \$\\partial_t u(t,x) + \\partial_x ( a(x) u(t,x) ) = 0\$
with boundary conditions `left_bc(t)`, `right_bc(t)`.

`D` is an SBP derivative operator, `Di` an associated dissipation operator
or `nothing`, `a(x)` the variable coefficient, and `split_form::Union{Val(false), Val(true)}`
determines whether the canonical split form or the conservative form should be
used.
"""
@auto_hash_equals struct VariableLinearAdvectionNonperiodicSemidiscretization{T, Derivative<:AbstractDerivativeOperator{T},
                                                                              Dissipation,
                                                                              SplitForm<:Union{Val{false}, Val{true}},
                                                                              LeftBC, RightBC} <: AbstractSemidiscretization
    derivative::Derivative
    dissipation::Dissipation
    a::Vector{T}
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm
    left_bc::LeftBC
    right_bc::RightBC

    function VariableLinearAdvectionNonperiodicSemidiscretization(derivative::Derivative, dissipation::Dissipation, afunc, split_form::SplitForm, left_bc::LeftBC, right_bc::RightBC) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}, LeftBC, RightBC}
        if dissipation !== nothing
            @argcheck size(derivative) == size(dissipation) DimensionMismatch
            @argcheck grid(derivative) == grid(dissipation) ArgumentError
        end
        N = size(derivative, 2)
        a = compute_coefficients(afunc, derivative)
        tmp1 = Array{T}(undef, N)
        tmp2 = Array{T}(undef, N)
        new{T,Derivative,Dissipation,SplitForm,LeftBC,RightBC}(derivative, dissipation, a, tmp1, tmp2, split_form, left_bc, right_bc)
    end
end


function Base.show(io::IO, semi::VariableLinearAdvectionNonperiodicSemidiscretization)
    if get(io, :compact, false)
        print(io, "Semidiscretization of the linear advection equation (nonperiodic)")
    else
        print(io, "Semidiscretization of the linear advection equation\n")
        print(io, "  \$ \\partial_t u(t,x) + \\partial_x ( a(x) u(t,x) ) = 0 \$ \n")
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


function godunov_flux_variablelinearadvection(uₗ::T, uᵣ::T, a::T) where {T<:Real}
    ifelse(a > 0, a*uₗ, a*uᵣ)
end

function (disc::VariableLinearAdvectionNonperiodicSemidiscretization)(du, u, p, t)
    @unpack a, tmp1, tmp2, derivative, dissipation, split_form, left_bc, right_bc = disc
    @boundscheck begin
        @argcheck length(u) == length(a)
        @argcheck length(u) == length(tmp1)
        @argcheck length(u) == length(tmp2)
        @argcheck length(u) == length(du)
    end

    # volume terms
    if typeof(split_form) <: Val{true}
        m1_2 = -one(eltype(u)) / 2

        ## a * D * u
        mul!(tmp1, derivative, u)
        @. du = m1_2 * a * tmp1
        ## u * D * a
        mul!(tmp1, derivative, a)
        @. du += m1_2 * u * tmp1
        ## D * a*u
        @. tmp2 = a*u
        mul!(tmp1, derivative, tmp2)
        @. du += m1_2 * tmp1
    else
        @. tmp2 = a * u
        mul!(tmp1, derivative, tmp2)
        @. du = -tmp1
    end

    # dissipation
    if dissipation !== nothing
        mul!(tmp1, dissipation, u)
        @. du += tmp1
    end

    # boundary conditions via Godunov's flux
    @inbounds fnum_left = godunov_flux_variablelinearadvection(left_bc(t), u[1], a[1])
    @inbounds du[1] += (fnum_left - a[1]*u[1]) / left_boundary_weight(derivative)
    @inbounds fnum_right = godunov_flux_variablelinearadvection(u[end], right_bc(t), a[end])
    @inbounds du[end] -= (fnum_right - a[end]*u[end]) / right_boundary_weight(derivative)

    nothing
end


"""
    VariableLinearAdvectionPeriodicSemidiscretization(D, Di, a, split_form)

A semidiscretization of the linear advection equation
    \$\\partial_t u(t,x) + \\partial_x ( a(x) u(t,x) ) = 0\$
with periodic boundary conditions.

`D` is a periodic SBP derivative operator, `Di` an associated dissipation operator
or `nothing`, `a(x)` the variable coefficient, and `split_form::Union{Val(false), Val(true)}`
determines whether the canonical split form or the conservative form should be
used.
"""
@auto_hash_equals struct VariableLinearAdvectionPeriodicSemidiscretization{T, Derivative<:AbstractPeriodicDerivativeOperator{T},
                                                                           Dissipation,
                                                                           SplitForm<:Union{Val{false}, Val{true}}} <: AbstractSemidiscretization
    derivative::Derivative
    dissipation::Dissipation
    a::Vector{T}
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm

    function VariableLinearAdvectionPeriodicSemidiscretization(derivative::Derivative, dissipation::Dissipation, afunc, split_form::SplitForm) where {T, Derivative<:AbstractPeriodicDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}}
        if dissipation !== nothing
            @argcheck size(derivative) == size(dissipation) DimensionMismatch
            @argcheck grid(derivative) == grid(dissipation) ArgumentError
        end
        N = size(derivative, 2)
        a = compute_coefficients(afunc, derivative)
        tmp1 = Array{T}(undef, N)
        tmp2 = Array{T}(undef, N)
        new{T,Derivative,Dissipation,SplitForm}(derivative, dissipation, a, tmp1, tmp2, split_form)
    end
end


function Base.show(io::IO, semi::VariableLinearAdvectionPeriodicSemidiscretization)
    if get(io, :compact, false)
        print(io, "Semidiscretization of the linear advection equation (periodic)")
    else
        print(io, "Semidiscretization of the linear advection equation\n")
        print(io, "  \$ \\partial_t u(t,x) + \\partial_x ( a(x) u(t,x) ) = 0 \$ \n")
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


function (disc::VariableLinearAdvectionPeriodicSemidiscretization)(du, u, p, t)
    @unpack a, tmp1, tmp2, derivative, dissipation, split_form = disc
    @boundscheck begin
        @argcheck length(u) == length(a)
        @argcheck length(u) == length(tmp1)
        @argcheck length(u) == length(tmp2)
        @argcheck length(u) == length(du)
    end

    # volume terms
    if typeof(split_form) <: Val{true}
        m1_2 = -one(eltype(u)) / 2

        ## a * D * u
        mul!(tmp1, derivative, u)
        @. du = m1_2 * a * tmp1
        ## u * D * a
        mul!(tmp1, derivative, a)
        @. du += m1_2 * u * tmp1
        ## D * a*u
        @. tmp2 = a*u
        mul!(tmp1, derivative, tmp2)
        @. du += m1_2 * tmp1
    else
        @. tmp2 = a * u
        mul!(tmp1, derivative, tmp2)
        @. du = -tmp1
    end

    # dissipation
    if dissipation !== nothing
        mul!(tmp1, dissipation, u)
        @. du += tmp1
    end

    nothing
end
