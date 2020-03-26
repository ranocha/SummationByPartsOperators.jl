
abstract type AbstractMesh1D end
abstract type AbstractPeriodicMesh1D <: AbstractMesh1D end

isperiodic(mesh::AbstractMesh1D) = false
isperiodic(mesh::AbstractPeriodicMesh1D) = true


function Base.show(io::IO, mesh::AbstractMesh1D)
  print(io,
    typeof(mesh), " with ", numcells(mesh), " cells in ", bounds(mesh))
end



"""
A uniform periodic mesh in one space dimension of `Nx` cells between
`xmin` and `xmax`.
"""
struct UniformPeriodicMesh1D{T<:Real} <: AbstractPeriodicMesh1D
  xmin::T
  xmax::T
  Nx::Int

  Δx::T

  function UniformPeriodicMesh1D{T}(xmin::T, xmax::T, Nx::Int) where T
    Nx > 0 || throw(ArgumentError("The number of elements `Nx` must be positive [`Nx == $Nx`]."))
    xmin < xmax || throw(ArgumentError("`xmin` must be smaller than `xmax` [`xmin == $xmin, xmax == $xmax`]."))

    new(xmin, xmax, Nx, (xmax-xmin)/Nx)
  end
end

function UniformPeriodicMesh1D(_xmin::Real, _xmax::Real, Nx::Integer)
  xmin, xmax = promote(_xmin, _xmax)
  UniformPeriodicMesh1D{typeof(xmin)}(xmin, xmax, Int(Nx))
end



"""
A uniform mesh in one space dimension of `Nx` cells between `xmin` and `xmax`.
"""
struct UniformMesh1D{T<:Real} <: AbstractMesh1D
  xmin::T
  xmax::T
  Nx::Int

  Δx::T

  function UniformMesh1D{T}(xmin::T, xmax::T, Nx::Int) where T
    Nx > 0 || throw(ArgumentError("The number of elements `Nx` must be positive [`Nx == $Nx`]."))
    xmin < xmax || throw(ArgumentError("`xmin` must be smaller than `xmax` [`xmin == $xmin, xmax == $xmax`]."))

    new(xmin, xmax, Nx, (xmax-xmin)/Nx)
  end
end

function UniformMesh1D(_xmin::Real, _xmax::Real, Nx::Integer)
  xmin, xmax = promote(_xmin, _xmax)
  UniformMesh1D{typeof(xmin)}(xmin, xmax, Int(Nx))
end

@inline numcells(mesh::Union{UniformMesh1D,UniformPeriodicMesh1D}) = mesh.Nx
@inline cell_indices(mesh::Union{UniformMesh1D,UniformPeriodicMesh1D}) = Base.OneTo(numcells(mesh))
@inline volume(cell::Int, mesh::Union{UniformMesh1D,UniformPeriodicMesh1D}) = mesh.Δx

function bounds(cell::Int, mesh::Union{UniformMesh1D,UniformPeriodicMesh1D})
  @unpack Δx = mesh

  xmin = mesh.xmin + (cell-1)*Δx
  xmin = nextfloat(xmin)

  xmax = mesh.xmin + cell*Δx
  xmax = prevfloat(xmax)

  xmin, xmax
end

@inline bounds(mesh::Union{UniformMesh1D,UniformPeriodicMesh1D}) = mesh.xmin, mesh.xmax



# struct CoupledOperator