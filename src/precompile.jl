
using SnoopPrecompile: @precompile_setup, @precompile_all_calls

@precompile_setup begin
  @static if !isdefined(Base, :get_extension)
    # This is required on Julia v1.8 since they there use Requires.jl to load
    # the code for StaticArrays.jl, see
    # https://github.com/ranocha/SummationByPartsOperators.jl/pull/165#issuecomment-1439550184

    # FIXME
    # This does not work:
    #   ERROR: LoadError: Evaluation into the closed module `StaticArrayInterface`
    #   breaks incremental compilation because the side effects will not be permanent.
    #   This is likely due to some other module mutating `StaticArrayInterface` with
    #   `eval` during precompilation - don't do this.
    # StaticArrayInterface.include(
    #   joinpath(dirname(dirname(pathof(StaticArrayInterface))),
    #            "ext", "StaticArrayInterfaceStaticArraysExt.jl"))

    # FIXME
    # This does not work:
    StaticArrayInterface.__init__()

    v = SVector(1, 2, 3)
    @info "While precompiling" static_length(v) typeof(static_length(v)) v
  end

  @precompile_all_calls begin
    let D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                         xmin=0.0, xmax=2.0, N=10)
      u = zeros(eltype(D), size(D, 2))
      D * u
    end

    let D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                                xmin=0.0, xmax=2.0, N=10)
      u = zeros(eltype(D), size(D, 2))
      D * u
    end
  end
end
