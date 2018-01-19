
doc"
    CubicPeriodicSemidiscretisation{T,Derivative,Dissipation}

A semidiscretisation of the cubic conservation law
    $\partial_t u(t,x) + \partial_x u(t,x)^3 = 0$
with periodic boundary conditions.
"
struct CubicPeriodicSemidiscretisation{T,Derivative<:AbstractDerivativeOperator{T},
                                        Dissipation<:AbstractDerivativeOperator{T},
                                        SplitForm<:Union{Val{false}, Val{true}}} <: AbstractSemidiscretisation
    derivative::Derivative
    dissipation::Dissipation
    tmp1::Vector{T}
    tmp2::Vector{T}
    split_form::SplitForm

    function CubicPeriodicSemidiscretisation(derivative::Derivative, dissipation::Dissipation, split_form::SplitForm=Val{false}()) where {T, Derivative<:AbstractDerivativeOperator{T}, Dissipation<:AbstractDerivativeOperator{T}, SplitForm<:Union{Val{false}, Val{true}}}
        @argcheck size(derivative) == size(dissipation) DimensionMismatch
        @argcheck grid(derivative) == grid(dissipation) ArgumentError
        N = size(derivative, 2)
        tmp1 = Array{T}(N)
        tmp2 = Array{T}(N)
        new{T,Derivative,Dissipation,SplitForm}(derivative, dissipation, tmp1, tmp2, split_form)
    end
end


function Base.show(io::IO, semidisc::CubicPeriodicSemidiscretisation)
    print(io, "Semidiscretisation of the cubic conservation law\n")
    print(io, "  \$ \partial_t u(t,x) + \partial_x u(t,x)^3 = 0 \$ \n")
    print(io, "using")
    if semidisc.split_form == Val{true}()
        print(io, " a split form and: \n")
    else
        print(io, " no split form and: \n")
    end
    print(io, semidisc.derivative)
    print(io, semidisc.dissipation)
end


function (disc::CubicPeriodicSemidiscretisation)(t, u, du)
    @unpack tmp1, tmp2, derivative, dissipation, split_form = disc
    @boundscheck begin
        @argcheck length(u) == length(tmp1)
        @argcheck length(u) == length(tmp2)
        @argcheck length(u) == length(du)
    end

    # volume terms
    if typeof(split_form) <: Val{true}
        mhalf = -one(eltype(u)) / 2

        ## u^2 * D * u
        @. tmp2 = u^2
        A_mul_B!(tmp1, derivative, u)
        @. du = mhalf * tmp2 * tmp1
        ## u * D * u^2
        A_mul_B!(tmp1, derivative, tmp2)
        @. du += mhalf * u * tmp1
        ## D * u^3
        @. tmp2 *= u
        A_mul_B!(tmp1, derivative, tmp2)
        @. du += mhalf * tmp1
    else
        @. tmp2 = u^3
        A_mul_B!(tmp1, derivative, tmp2)
        @. du = -tmp1
    end

    # dissipation
    A_mul_B!(tmp1, dissipation, u)
    @. du += tmp1

    nothing
end

