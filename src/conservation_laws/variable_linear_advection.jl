"""
    VariableLinearAdvectionNonperiodicSemidiscretisation

A semidiscretisation of the linear advection equation
    \$\\partial_t u(t,x) + \\partial_x ( a(x) u(t,x) ) = 0\$
with boundary conditions `left_bc(t)`, `right_bc(t)`.
"""
struct VariableLinearAdvectionNonperiodicSemidiscretisation{T,Derivative<:AbstractDerivativeOperator{T},
                                                            Dissipation,
                                                            SplitForm<:Union{Val{false}, Val{true}},
                                                            LeftBC, RightBC} <: AbstractSemidiscretisation
    derivative::Derivative
    dissipation::Dissipation
    a::Vector{T}
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm
    left_bc::LeftBC
    right_bc::RightBC

    function VariableLinearAdvectionNonperiodicSemidiscretisation(derivative::Derivative, dissipation::Dissipation, afunc, split_form::SplitForm, left_bc::LeftBC, right_bc::RightBC) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}, LeftBC, RightBC}
        if dissipation != nothing
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


function Base.show(io::IO, semidisc::VariableLinearAdvectionNonperiodicSemidiscretisation)
    print(io, "Semidiscretisation of the linear advection equation\n")
    print(io, "  \$ \\partial_t u(t,x) + \\partial_x ( a(x) u(t,x) ) = 0 \$ \n")
    print(io, "with nonperiodic boundaries using")
    if semidisc.split_form == Val{true}()
        print(io, " a split form and: \n")
    else
        print(io, " no split form and: \n")
    end
    print(io, semidisc.derivative)
    print(io, semidisc.dissipation)
end


function godunov_flux_variablelinearadvection(uₗ::T, uᵣ::T, a::T) where {T<:Real}
    ifelse(a > 0, a*uₗ, a*uᵣ)
end

function (disc::VariableLinearAdvectionNonperiodicSemidiscretisation)(du, u, p, t)
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
    if dissipation != nothing
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
