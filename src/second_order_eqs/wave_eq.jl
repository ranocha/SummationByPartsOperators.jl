"""
    WaveEquationNonperiodicSemidiscretization(D, left_bc, right_bc)

A semidiscretization of the linear wave equation
    \$\\partial_t^2 u(t,x) = \\partial_x^2 u(t,x)\$.

`D` is assumed to be a second-derivative SBP operator and the boundary conditions
can be `Val(:HomogeneousNeumann)`, `Val(:HomogeneousDirichlet)`, or
`Val(:NonReflecting)`.
"""
@auto_hash_equals struct WaveEquationNonperiodicSemidiscretization{Derivative<:AbstractDerivativeOperator,
                                                                   LeftBC, RightBC} <: AbstractSemidiscretization
    derivative::Derivative
    left_bc::LeftBC
    right_bc::RightBC
end

function semidiscretize(du0func, u0func, semi::WaveEquationNonperiodicSemidiscretization, tspan)
  du0 = compute_coefficients(du0func, semi.derivative)
  u0 = compute_coefficients(u0func, semi.derivative)
  iip = true # is in-place
  return SecondOrderODEProblem{iip}(semi, du0, u0, tspan)
end


function Base.show(io::IO, semi::WaveEquationNonperiodicSemidiscretization)
    if get(io, :compact, false)
        print(io, "Semidiscretization of the linear wave equation (non-periodic)")
    else
        print(io, "Semidiscretization of the linear wave equation\n")
        print(io, "  \$ \\partial_t^2 u(t,x) = \\partial_x^2 u(t,x) \$ \n")
        print(io, "with nonperiodic boundaries using")
        print(io, semi.derivative)
    end
end


function (disc::WaveEquationNonperiodicSemidiscretization)(ddu, du, u, p, t)
    @unpack derivative, left_bc, right_bc = disc
    @boundscheck begin
        @argcheck length(u) == length(du)
        @argcheck length(u) == length(ddu)
    end

    # volume terms of the second derivative
    mul!(ddu, derivative, u)

    # boundary conditions using SATs
    if left_bc == Val(:HomogeneousNeumann)
        @inbounds ddu[1] += derivative_left(derivative, u, Val(1)) / left_boundary_weight(derivative)
    elseif left_bc == Val(:HomogeneousDirichlet)
        mul_transpose_derivative_left!(ddu, derivative, Val(1), -u[1] / left_boundary_weight(derivative), true)
        # only a rough estimate, should probably be improved
        # see Mattsson, Ham, Iaccarino (2009) Stable boundary ...
        # parameter α in Table 1 instead of accuracy_order(derivative)
        @inbounds ddu[1] -= accuracy_order(derivative) * u[1] / left_boundary_weight(derivative)^2
    elseif left_bc == Val(:NonReflecting)
        @inbounds ddu[1] -= (du[1] - derivative_left(derivative, u, Val(1))) / left_boundary_weight(derivative)
    else
        throw(ArgumentError("Boundary condition $left_bc not implemented."))
    end

    if right_bc == Val(:HomogeneousNeumann)
       @inbounds ddu[end] -= derivative_right(derivative, u, Val(1)) / right_boundary_weight(derivative)
    elseif right_bc == Val(:HomogeneousDirichlet)
        mul_transpose_derivative_right!(ddu, derivative, Val(1), u[end] / right_boundary_weight(derivative), true)
        # only a rough estimate, should probably be improved, see above
        @inbounds ddu[end] -= accuracy_order(derivative) * u[end] / right_boundary_weight(derivative)^2
    elseif right_bc == Val(:NonReflecting)
        @inbounds ddu[end] -= (du[end] + derivative_right(derivative, u, Val(1))) / right_boundary_weight(derivative)
    else
        throw(ArgumentError("Boundary condition $right_bc not implemented."))
    end

    nothing
end
