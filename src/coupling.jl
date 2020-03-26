
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



struct UniformMeshGrid1D{T, Mesh<:Union{UniformMesh1D,UniformPeriodicMesh1D}, Grid<:AbstractVector{T}} <: AbstractArray{T,1}
  mesh::Mesh
  grid::Grid
  continuous::Bool
end

isperiodic(meshgrid::UniformMeshGrid1D) = isperiodic(meshgrid.mesh)
numcells(meshgrid::UniformMeshGrid1D) = numcells(meshgrid.mesh)
bounds(meshgrid::UniformMeshGrid1D) = bounds(meshgrid.mesh)
bounds(i::Int, meshgrid::UniformMeshGrid1D) = bounds(i, meshgrid.mesh)
iscontinuous(meshgrid::UniformMeshGrid1D) = meshgrid.continuous

function Base.show(io::IO, meshgrid::UniformMeshGrid1D)
  print(io,
    "UniformMeshGrid1D with ", numcells(meshgrid), " cells with ", length(meshgrid.grid), " nodes in ", bounds(meshgrid))
end

function Base.length(meshgrid::UniformMeshGrid1D)
  if iscontinuous(meshgrid)
    if isperiodic(meshgrid)
      numcells(meshgrid) * (length(meshgrid.grid) - 1)
    else
      numcells(meshgrid) * (length(meshgrid.grid) - 1) + 1
    end
  else
    numcells(meshgrid) * length(meshgrid.grid)
  end
end

function Base.getindex(meshgrid::UniformMeshGrid1D, i::Int)
  N = length(meshgrid)
  @boundscheck begin
    @argcheck i > 0
    @argcheck i <= N
  end
  if iscontinuous(meshgrid)
    num_nodes_per_cell = length(meshgrid.grid) - 1
    if i <= num_nodes_per_cell + 1
      cell = 1
      node = i
    else
      cell = (i - 2) ÷ num_nodes_per_cell + 1
      node = i - (cell-1) * num_nodes_per_cell
    end
  else
    num_nodes_per_cell = length(meshgrid.grid)
    cell = (i-1) ÷ num_nodes_per_cell + 1
    node = i - (cell-1) * num_nodes_per_cell
  end
  xmin, xmax = bounds(cell, meshgrid)
  ymin, ymax = first(meshgrid.grid), last(meshgrid.grid)
  (meshgrid.grid[node] - ymin) * (xmax - xmin) / (ymax - ymin) + xmin
end

Base.size(meshgrid::UniformMeshGrid1D) = (length(meshgrid),)


struct UniformCoupledOperator{T, Dtype<:AbstractNonperiodicDerivativeOperator{T}, MeshGrid <: UniformMeshGrid1D{T}} <: AbstractNonperiodicDerivativeOperator{T}
  D::Dtype
  meshgrid::MeshGrid

  function UniformCoupledOperator(D::Dtype, mesh::Mesh, continuous::Bool) where {T, Dtype<:AbstractNonperiodicDerivativeOperator{T}, Mesh<:AbstractMesh1D}
    meshgrid = UniformMeshGrid1D(mesh, grid(D), continuous)
    new{T, Dtype, typeof(meshgrid)}(D, meshgrid)
  end
end
