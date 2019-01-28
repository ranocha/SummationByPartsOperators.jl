"""
WaveEquationNonperiodicSemidiscretisation

A semidiscretisation of the linear wave equation
    \$\\partial_t^2 u(t,x) = \\partial_x^2 u(t,x)\$.
"""
struct WaveEquationNonperiodicSemidiscretisation{Derivative<:AbstractDerivativeOperator,
                                                 LeftBC, RightBC} <: AbstractSemidiscretisation
    derivative::Derivative
    left_bc::LeftBC
    right_bc::RightBC
end

function semidiscretise(du0func, u0func, semidisc::WaveEquationNonperiodicSemidiscretisation, tspan)
  du0 = compute_coefficients(du0func, semidisc.derivative)
  u0 = compute_coefficients(u0func, semidisc.derivative)
  ode = SecondOrderODEProblem(semidisc, du0, u0, tspan)
end


function Base.show(io::IO, semidisc::WaveEquationNonperiodicSemidiscretisation)
    print(io, "Semidiscretisation of the linear wave equation\n")
    print(io, "  \$ \\partial_t^2 u(t,x) = \\partial_x^2 u(t,x) \$ \n")
    print(io, "with nonperiodic boundaries using")
    print(io, semidisc.derivative)
end


function (disc::WaveEquationNonperiodicSemidiscretisation)(ddu, du, u, p, t)
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
        add_transpose_derivative_left!(ddu, derivative, Val(1), -u[1] / left_boundary_weight(derivative))
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
        add_transpose_derivative_right!(ddu, derivative, Val(1), u[end] / right_boundary_weight(derivative))
        # only a rough estimate, should probably be improved, see above
        @inbounds ddu[end] -= accuracy_order(derivative) * u[end] / right_boundary_weight(derivative)^2
    elseif right_bc == Val(:NonReflecting)
        @inbounds ddu[end] -= (du[end] + derivative_right(derivative, u, Val(1))) / right_boundary_weight(derivative)
    else
        throw(ArgumentError("Boundary condition $right_bc not implemented."))
    end

    nothing
end
