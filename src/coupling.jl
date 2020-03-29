
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

function left_cell(cell::Int, mesh::AbstractMesh1D)
  N = numcells(mesh)
  @boundscheck begin
    @assert(1 <= cell <= N)
  end
  isperiodic(mesh) && cell == 1 ? N : cell-1
end

function right_cell(cell::Int, mesh::AbstractMesh1D)
  N = numcells(mesh)
  @boundscheck begin
    @assert(1 <= cell <= N)
  end
  isperiodic(mesh) && cell == N ? 1 : cell+1
end



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
    "UniformMeshGrid1D with ", numcells(meshgrid), " cells with ", length(meshgrid.grid), " nodes in ", bounds(meshgrid), "\n ")
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



struct UniformCoupledOperator{T, Dtype<:AbstractNonperiodicDerivativeOperator{T}, MeshGrid <: UniformMeshGrid1D{T}, Coupling<:Union{Val{:continuous}, Val{:plus}, Val{:central}, Val{:minus}}} <: AbstractNonperiodicDerivativeOperator{T}
  D::Dtype
  meshgrid::MeshGrid
  coupling::Coupling

  function UniformCoupledOperator(D::Dtype, mesh::Mesh, coupling::Coupling) where {T, Dtype<:AbstractNonperiodicDerivativeOperator{T}, Mesh<:AbstractMesh1D, Coupling<:Union{Val{:continuous}, Val{:plus}, Val{:central}, Val{:minus}}}
    meshgrid = UniformMeshGrid1D(mesh, grid(D), coupling===Val(:continuous))
    if derivative_order(D) != 1
      throw(ArgumentError("Not implemented yet"))
    end
    new{T, Dtype, typeof(meshgrid), Coupling}(D, meshgrid, coupling)
  end
end

iscontinuous(cD::UniformCoupledOperator) = iscontinuous(cD.meshgrid)

function Base.show(io::IO, cD::UniformCoupledOperator)
  print(io, cD.D)
  if iscontinuous(cD)
    print(io, "coupled continuously")
  else
    print(io, "coupled discontinuously (upwind: ", cD.coupling, ")")
  end
  print(io, " on the mesh \n")
  print(io, cD.meshgrid.mesh, "\n ")
end

Base.eltype(cD::UniformCoupledOperator) = eltype(cD.D)
derivative_order(cD::UniformCoupledOperator) = derivative_order(cD.D)
accuracy_order(cD::UniformCoupledOperator) = accuracy_order(cD.D)
source_of_coefficients(cD::UniformCoupledOperator) = source_of_coefficients(cD.D)

function left_boundary_weight(cD::UniformCoupledOperator)
  @unpack D, meshgrid = cD
  @unpack mesh, grid = meshgrid
  ymin, ymax = first(grid), last(grid)
  cell = 1
  xmin, xmax = bounds(cell, mesh)
  factor = (xmax - xmin) / (ymax - ymin)
  factor * left_boundary_weight(D)
end

function right_boundary_weight(cD::UniformCoupledOperator)
  @unpack D, meshgrid = cD
  @unpack mesh, grid = meshgrid
  ymin, ymax = first(grid), last(grid)
  cell = numcells(mesh)
  xmin, xmax = bounds(cell, mesh)
  factor = (xmax - xmin) / (ymax - ymin)
  factor * right_boundary_weight(D)
end

function couple_continuosly(D::AbstractNonperiodicDerivativeOperator, mesh::AbstractMesh1D)
  UniformCoupledOperator(D, mesh, Val(:continuous))
end

function couple_discontinuosly(D::AbstractNonperiodicDerivativeOperator, mesh::AbstractMesh1D, coupling::Union{Val{:plus}, Val{:central}, Val{:minus}}=Val(:central))
  UniformCoupledOperator(D, mesh, coupling)
end

grid(cD::UniformCoupledOperator) = cD.meshgrid


function mass_matrix(cD::UniformCoupledOperator)
  m = ones(eltype(cD), size(cD, 1))
  scale_by_mass_matrix!(m, cD)
  Diagonal(vec(m))
end

for (fname, op) in ((:scale_by_mass_matrix!, Base.:*), (:scale_by_inverse_mass_matrix!, Base.:/))
  @eval function $fname(u::AbstractVector, cD::UniformCoupledOperator, factor=true)
    @unpack D, meshgrid = cD
    @unpack mesh, grid = meshgrid
    ymin, ymax = first(grid), last(grid)

    if iscontinuous(meshgrid)
      num_nodes_per_cell = length(grid) - 1

      cell = 1
      xmin, xmax = bounds(cell, mesh)
      jac = (xmax - xmin) / (ymax - ymin)
      if isperiodic(mesh)
        next_cell = left_cell(cell, mesh)
        xmin, xmax = bounds(next_cell, mesh)
        next_jac = (xmax - xmin) / (ymax - ymin)
        u[1] *= $op(factor, jac * get_weight(D, 1) + next_jac * get_weight(D, num_nodes_per_cell+1))
      else
        u[1] *= $op(factor, jac * get_weight(D, 1))
      end
      for i in 2:num_nodes_per_cell
        u[i] *= $op(factor, jac * get_weight(D, i))
      end
      if !isperiodic(mesh) && numcells(mesh) == 1
        u[num_nodes_per_cell+1] *= $op(factor, jac * get_weight(D, num_nodes_per_cell+1))
      end
      if numcells(mesh) == 1
        return u
      end

      for cell in 2:numcells(mesh)-1
        xmin, xmax = bounds(cell, mesh)
        jac = (xmax - xmin) / (ymax - ymin)
        next_cell = left_cell(cell, mesh)
        xmin, xmax = bounds(next_cell, mesh)
        next_jac = (xmax - xmin) / (ymax - ymin)
        u[(cell-1)*num_nodes_per_cell+1] *= $op(factor, jac * get_weight(D, 1) + next_jac * get_weight(D, num_nodes_per_cell+1))
        for i in 2:num_nodes_per_cell
          u[(cell-1)*num_nodes_per_cell+i] *= $op(factor, jac * get_weight(D, i))
        end
      end

      cell = numcells(mesh)
      xmin, xmax = bounds(cell, mesh)
      jac = (xmax - xmin) / (ymax - ymin)
      next_cell = left_cell(cell, mesh)
      xmin, xmax = bounds(next_cell, mesh)
      next_jac = (xmax - xmin) / (ymax - ymin)
      u[(cell-1)*num_nodes_per_cell+1] *= $op(factor, jac * get_weight(D, 1) + next_jac * get_weight(D, num_nodes_per_cell+1))
      for i in 2:num_nodes_per_cell
        u[(cell-1)*num_nodes_per_cell+i] *= $op(factor, jac * get_weight(D, i))
      end
      if !isperiodic(mesh)
        u[(cell-1)*num_nodes_per_cell+num_nodes_per_cell+1] *= $op(factor, jac * get_weight(D, num_nodes_per_cell+1))
      end
    else # discontinuous coupling
      num_nodes_per_cell = length(grid)
      for cell in 1:numcells(mesh)
        xmin, xmax = bounds(cell, mesh)
        jac = (xmax - xmin) / (ymax - ymin)
        idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell
        $fname(view(u, idx), D, $op(factor, jac))
      end
    end

    u
  end
end

function integrate(func, u::AbstractVector, cD::UniformCoupledOperator)
  @unpack D, meshgrid = cD
  @unpack mesh, grid = meshgrid
  ymin, ymax = first(grid), last(grid)

  if iscontinuous(meshgrid)
    num_nodes_per_cell = length(grid) - 1

    cell = 1
    xmin, xmax = bounds(cell, mesh)
    jac = (xmax - xmin) / (ymax - ymin)
    if isperiodic(mesh)
      next_cell = left_cell(cell, mesh)
      xmin, xmax = bounds(next_cell, mesh)
      next_jac = (xmax - xmin) / (ymax - ymin)
      res = func(u[1]) * (jac * get_weight(D, 1) + next_jac * get_weight(D, num_nodes_per_cell+1))
    else
      res = func(u[1]) * (jac * get_weight(D, 1))
    end
    for i in 2:num_nodes_per_cell
      res += func(u[i]) * (jac * get_weight(D, i))
    end
    if !isperiodic(mesh) && numcells(mesh) == 1
      res += func(u[num_nodes_per_cell+1]) * (jac * get_weight(D, num_nodes_per_cell+1))
    end
    if numcells(mesh) == 1
      return res
    end

    for cell in 2:numcells(mesh)-1
      xmin, xmax = bounds(cell, mesh)
      jac = (xmax - xmin) / (ymax - ymin)
      next_cell = left_cell(cell, mesh)
      xmin, xmax = bounds(next_cell, mesh)
      next_jac = (xmax - xmin) / (ymax - ymin)
      res += func(u[(cell-1)*num_nodes_per_cell+1]) * (jac * get_weight(D, 1) + next_jac * get_weight(D, num_nodes_per_cell+1))
      for i in 2:num_nodes_per_cell
        res += func(u[(cell-1)*num_nodes_per_cell+i]) * (jac * get_weight(D, i))
      end
    end

    cell = numcells(mesh)
    xmin, xmax = bounds(cell, mesh)
    jac = (xmax - xmin) / (ymax - ymin)
    next_cell = left_cell(cell, mesh)
    xmin, xmax = bounds(next_cell, mesh)
    next_jac = (xmax - xmin) / (ymax - ymin)
    res += func(u[(cell-1)*num_nodes_per_cell+1]) * (jac * get_weight(D, 1) + next_jac * get_weight(D, num_nodes_per_cell+1))
    for i in 2:num_nodes_per_cell
      res += func(u[(cell-1)*num_nodes_per_cell+i]) * (jac * get_weight(D, i))
    end
    if !isperiodic(mesh)
      res += func(u[(cell-1)*num_nodes_per_cell+num_nodes_per_cell+1]) * (jac * get_weight(D, num_nodes_per_cell+1))
    end
  else # discontinuous coupling
    num_nodes_per_cell = length(grid)
    cell = 1
    xmin, xmax = bounds(cell, mesh)
    jac = (xmax - xmin) / (ymax - ymin)
    idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell
    res = jac * integrate(view(u, idx), D)
    for cell in 2:numcells(mesh)
      xmin, xmax = bounds(cell, mesh)
      jac = (xmax - xmin) / (ymax - ymin)
      idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell
      res += jac * integrate(view(u, idx), D)
    end
  end

  res
end


function mul!(dest::AbstractVector, cD::UniformCoupledOperator, u::AbstractVector, α=true)
  N, _ = size(cD)
  @boundscheck begin
    @argcheck N == length(u)
    @argcheck N == length(dest)
  end

  @unpack D, meshgrid, coupling = cD
  if coupling === Val(:continuous)
    mul!(dest, D, meshgrid, coupling, u, α)
    scale_by_inverse_mass_matrix!(dest, cD)
  else
    mul!(dest, D, meshgrid, coupling, u, α)
  end
  dest
end

function mul!(_dest::AbstractVector, D::AbstractNonperiodicDerivativeOperator, meshgrid::UniformMeshGrid1D, coupling::Union{Val{:plus}, Val{:central}, Val{:minus}}, _u::AbstractVector, α=true)
  @unpack mesh, grid = meshgrid
  dest = reshape(_dest, length(grid), numcells(mesh))
  u    = reshape(_u,    length(grid), numcells(mesh))
  ymin, ymax = first(grid), last(grid)
  half = one(eltype(D)) / 2

  cell = 1
  xmin, xmax = bounds(cell, mesh)
  factor = (ymax - ymin) / (xmax - xmin)
  mul!(view(dest, :, cell), D, view(u, :, cell), α*factor)
  if numcells(mesh) == 1 && !isperiodic(mesh)
    return _dest
  end
  if coupling === Val(:central)
    dest[end, cell] += half * (u[1, right_cell(cell, mesh)] - u[end, cell]) * α*factor / right_boundary_weight(D)
  elseif coupling === Val(:plus)
    dest[end, cell] +=        (u[1, right_cell(cell, mesh)] - u[end, cell]) * α*factor / right_boundary_weight(D)
  end
  if isperiodic(mesh)
    if coupling === Val(:central)
      dest[1  , cell] += half * (u[1, cell] - u[end,  left_cell(cell, mesh)]) * α*factor / left_boundary_weight(D)
    elseif coupling === Val(:minus)
      dest[1  , cell] +=        (u[1, cell] - u[end,  left_cell(cell, mesh)]) * α*factor / left_boundary_weight(D)
    end
  end
  if numcells(mesh) == 1
    return _dest
  end

  for cell in 2:numcells(mesh)-1
    xmin, xmax = bounds(cell, mesh)
    factor = (ymax - ymin) / (xmax - xmin)
    mul!(view(dest, :, cell), D, view(u, :, cell), α*factor)
    if coupling === Val(:central)
      dest[1  , cell] += half * (u[1, cell] - u[end,  left_cell(cell, mesh)]) * α*factor / left_boundary_weight(D)
      dest[end, cell] += half * (u[1, right_cell(cell, mesh)] - u[end, cell]) * α*factor / right_boundary_weight(D)
    elseif coupling === Val(:plus)
      dest[end, cell] +=        (u[1, right_cell(cell, mesh)] - u[end, cell]) * α*factor / right_boundary_weight(D)
    elseif coupling === Val(:minus)
      dest[1  , cell] +=        (u[1, cell] - u[end,  left_cell(cell, mesh)]) * α*factor / left_boundary_weight(D)
    end
  end

  cell = numcells(mesh)
  xmin, xmax = bounds(cell, mesh)
  factor = (ymax - ymin) / (xmax - xmin)
  mul!(view(dest, :, cell), D, view(u, :, cell), α*factor)
  if coupling === Val(:central)
    dest[1  , cell] += half * (u[1, cell] - u[end,  left_cell(cell, mesh)]) * α*factor / left_boundary_weight(D)
  elseif coupling === Val(:minus)
    dest[1  , cell] +=        (u[1, cell] - u[end,  left_cell(cell, mesh)]) * α*factor / left_boundary_weight(D)
  end
  if isperiodic(mesh)
    if coupling === Val(:central)
      dest[end, cell] += half * (u[1, right_cell(cell, mesh)] - u[end, cell]) * α*factor / right_boundary_weight(D)
    elseif coupling === Val(:plus)
      dest[end, cell] +=        (u[1, right_cell(cell, mesh)] - u[end, cell]) * α*factor / right_boundary_weight(D)
    end
  end

  _dest
end

function mul!(dest::AbstractVector, D::AbstractNonperiodicDerivativeOperator, meshgrid::UniformMeshGrid1D, coupling::Val{:continuous}, u::AbstractVector, α=true)
  @unpack mesh, grid = meshgrid
  ymin, ymax = first(grid), last(grid)
  num_nodes_per_cell = length(grid) - 1
  # TODO: remove these allocations?
  utmp = similar(u, length(grid))
  desttmp = similar(dest, length(grid))

  cell = 1
  xmin, xmax = bounds(cell, mesh)
  factor = (ymax - ymin) / (xmax - xmin)
  if isperiodic(mesh) && numcells(mesh) == 1
    utmp[1:end-1] .= u
    utmp[end]      = u[1]
    mul!(desttmp, D, utmp, α*factor)
    scale_by_mass_matrix!(desttmp, D, inv(factor))
    @. dest = desttmp[1:end-1]
    dest[1] += desttmp[end]
  else
    idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell+1
    utmp .= u[idx]
    mul!(desttmp, D, utmp, α*factor)
    scale_by_mass_matrix!(desttmp, D, inv(factor))
    @. dest[idx] = desttmp
  end
  if numcells(mesh) == 1
    return dest
  end

  for cell in 2:numcells(mesh)-1
    xmin, xmax = bounds(cell, mesh)
    factor = (ymax - ymin) / (xmax - xmin)
    idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell+1
    utmp .= u[idx]
    mul!(desttmp, D, utmp, α*factor)
    scale_by_mass_matrix!(desttmp, D, inv(factor))
    dest[idx[1]] += desttmp[1]
    @. dest[idx[2:end]] = desttmp[2:end]
  end

  cell = numcells(mesh)
  xmin, xmax = bounds(cell, mesh)
  factor = (ymax - ymin) / (xmax - xmin)
  if isperiodic(mesh)
    idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell
    utmp[1:end-1] .= u[idx]
    utmp[end]      = u[1]
    mul!(desttmp, D, utmp, α*factor)
    scale_by_mass_matrix!(desttmp, D, inv(factor))
    dest[idx[1]] += desttmp[1]
    @. dest[idx[2:end]] = desttmp[2:end-1]
    dest[1] += desttmp[end]
  else
    idx = (cell-1)*num_nodes_per_cell+1:cell*num_nodes_per_cell+1
    utmp .= u[idx]
    mul!(desttmp, D, utmp, α*factor)
    scale_by_mass_matrix!(desttmp, D, inv(factor))
    dest[idx[1]] += desttmp[1]
    @. dest[idx[2:end]] = desttmp[2:end]
  end

  dest
end
