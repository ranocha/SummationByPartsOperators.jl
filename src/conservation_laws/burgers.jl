
doc"
BurgersPeriodicSemidiscretisation{T,Derivative,Dissipation}

A semidiscretisation of Burgers' equation
$\partial_t u(t,x) + \partial_x \frac{u(t,x)^2}{2} = 0$
with periodic boundary conditions.
"
struct BurgersPeriodicSemidiscretisation{T,Derivative<:AbstractDerivativeOperator{T},
                                            Dissipation<:AbstractDerivativeOperator{T},
                                            SplitForm<:Union{Val{false}, Val{true}}} <: AbstractSemidiscretisation
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm

    function BurgersPeriodicSemidiscretisation(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm=Val{false}()) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation<:AbstractDerivativeOperator{T}, SplitForm<:Union{Val{false}, Val{true}}}
        @argcheck size(derivative) == size(dissipation) DimensionMismatch
        @argcheck grid(derivative) == grid(dissipation) ArgumentError
        N = size(derivative, 2)
        tmp1 = Array{T}(N)
        tmp2 = Array{T}(N)
        new{T,Derivative,Dissipation,SplitForm}(derivative, dissipation, tmp1, tmp2, split_form)
    end
end


function Base.show(io::IO, semidisc::BurgersPeriodicSemidiscretisation)
    print(io, "Semidiscretisation of Burgers' equation\n")
    print(io, "  \$ \partial_t u(t,x) + \partial_x \frac{u(t,x)^2}{2} = 0 \$ \n")
    print(io, "using")
    if semidisc.split_form == Val{true}()
        print(io, " a split form and: \n")
    else
        print(io, " no split form and: \n")
    end
    print(io, semidisc.derivative)
    print(io, semidisc.dissipation)
end


function (disc::BurgersPeriodicSemidiscretisation)(t, u, du)
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
        A_mul_B!(tmp1, derivative, u)
        @. du = m1_3 * u * tmp1
        ## D * u^2
        @. tmp2 = u^2
        A_mul_B!(tmp1, derivative, tmp2)
        @. du += m1_3 * tmp1
    else
        m1_2 = -one(eltype(u)) / 2
        @. tmp2 = u^2
        A_mul_B!(tmp1, derivative, tmp2)
        @. du = m1_2 * tmp1
    end

    # dissipation
    A_mul_B!(tmp1, dissipation, u)
    @. du += tmp1

    nothing
end

