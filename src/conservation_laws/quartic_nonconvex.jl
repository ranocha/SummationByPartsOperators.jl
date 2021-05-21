
"""
    QuarticNonconvexPeriodicSemidiscretisation(D, Di, split_form)

A semidiscretisation of the quartic nonconvex conservation law
    \$\\partial_t u(t,x) + \\partial_x ( u(t,x)^4 - 10 u(t,x)^2 + 3 u(t,x) ) = 0\$
with periodic boundary conditions.

`D` is a first-derivative SBP operator, `Di` an associated dissipation operator
or `nothing`, and `split_form::Union{Val(true), Val(false)}` determines whether
the canonical split form or the conservative form is used.
"""
struct QuarticNonconvexPeriodicSemidiscretisation{T,Derivative<:AbstractDerivativeOperator{T},
                                                  Dissipation,
                                                  SplitForm<:Union{Val{false}, Val{true}}} <: AbstractSemidiscretisation
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm

    function QuarticNonconvexPeriodicSemidiscretisation(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm=Val{false}()) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation, SplitForm<:Union{Val{false}, Val{true}}}
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


function Base.show(io::IO, semidisc::QuarticNonconvexPeriodicSemidiscretisation)
    print(io, "Semidiscretisation of the quartic nonconvex conservation law\n")
    print(io, "  \$ \\partial_t u(t,x) + \\partial_x ( u(t,x)^4 - 10 u(t,x)^2 + 3 u(t,x) ) = 0 \$ \n")
    print(io, "with periodic boundaries using")
    if semidisc.split_form == Val{true}()
        print(io, " a split form and: \n")
    else
        print(io, " no split form and: \n")
    end
    print(io, semidisc.derivative)
    print(io, semidisc.dissipation)
end


function (disc::QuarticNonconvexPeriodicSemidiscretisation)(du, u, p, t)
    @unpack tmp1, tmp2, derivative, dissipation, split_form = disc
    @boundscheck begin
        @argcheck length(u) == length(tmp1)
        @argcheck length(u) == length(tmp2)
        @argcheck length(u) == length(du)
    end

    # volume terms
    if typeof(split_form) <: Val{true}
        m_2_5 = -2*one(eltype(u)) / 5
        m_20_3 = -20*one(eltype(u)) / 3
        m_3 = -3*one(eltype(u))

        # # D u^4
        @. tmp2 = u^4
        mul!(tmp1, derivative, tmp2)
        @. du = m_2_5 * tmp1

        @. tmp2 = u^3
        mul!(tmp1, derivative, tmp2)
        @. du += m_2_5 * u * tmp1

        @. tmp2 = u^2
        mul!(tmp1, derivative, tmp2)
        @. du += m_2_5 * u^2 * tmp1

        mul!(tmp1, derivative, u)
        @. du += m_2_5 * u^3 * tmp1

        # D u^2
        @. tmp2 = u^2
        mul!(tmp1, derivative, tmp2)
        @. du -= m_20_3 * tmp1

        mul!(tmp1, derivative, u)
        @. du -= m_20_3 * u * tmp1

        # D u
        mul!(tmp1, derivative, u)
        @. du += m_3 * tmp1
    else
        @. tmp2 = u^4 - 10 * u^2 + 3 * u
        mul!(tmp1, derivative, tmp2)
        @. du = -tmp1
    end

    # dissipation
    if dissipation != nothing
        mul!(tmp1, dissipation, u)
        @. du += tmp1
    end

    nothing
end
