using SnoopPrecompile: @precompile_all_calls

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
