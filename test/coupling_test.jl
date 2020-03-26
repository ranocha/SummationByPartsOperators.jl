using Test
using LinearAlgebra
using SummationByPartsOperators

@testset "Coupled Legendre operators" begin
  for T in (Float32, Float64)
    for degree in 1:5
      ymin = T(0)
      ymax = T(2)
      D = legendre_derivative_operator(ymin, ymax, degree+1)

      for N in 1:3
        xmin = T(-1)
        xmax = T(2)
        mesh = UniformMesh1D(xmin, xmax, N)
        cD_continuous = couple_continuosly(D, mesh)
        cD_central    = couple_discontinuosly(D, mesh)
        cD_plus       = couple_discontinuosly(D, mesh, Val(:plus))
        cD_minus      = couple_discontinuosly(D, mesh, Val(:minus))
        for cD in (cD_continuous, cD_central, cD_plus, cD_minus)
          print(devnull, cD)
          x = grid(cD)
          @test norm(cD * x.^0) < 100N * eps(T)
          for k in 1:degree
            @test cD * x.^k ≈ k .* x.^(k-1)
          end

          u = sinpi.(x)
          v = copy(u)
          SummationByPartsOperators.scale_by_mass_matrix!(u, cD)
          SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, cD)
          @test u ≈ v
        end
        M = mass_matrix(cD_central)
        @test M ≈ mass_matrix(cD_plus)
        @test M ≈ mass_matrix(cD_minus)
        @test sum(M) ≈ xmax - xmin
        Dp = Matrix(cD_plus)
        Dm = Matrix(cD_minus)
        diss = M * (Dp - Dm)
        @test diss ≈ diss'
        @test maximum(eigvals(Symmetric(diss))) < 100N * eps(T)
        Dc = Matrix(cD_central)
        res = M * Dc + Dc' * M
        res[1, 1] += 1
        res[end, end] -= 1
        @test norm(res) < 100N * eps(T)
        res = M * Dp + Dm' * M
        res[1, 1] += 1
        res[end, end] -= 1
        @test norm(res) < 100N * eps(T)

        cD1 = couple_discontinuosly(cD_central, mesh)
        cD2 = couple_discontinuosly(D, UniformMesh1D(xmin, xmax, N^2))
        @test grid(cD1) ≈ grid(cD2)
        @test mass_matrix(cD1) ≈ mass_matrix(cD2)
        @test Matrix(cD1) ≈ Matrix(cD2)

        Mcont = mass_matrix(cD_continuous)
        @test sum(Mcont) ≈ xmax - xmin
        Dcont = Matrix(cD_continuous)
        res = Mcont * Dcont + Dcont' * Mcont
        res[1, 1] += 1
        res[end, end] -= 1
        @test norm(res) < 100N * eps(T)
      end

      for N in 1:3
        xmin = T(-1)
        xmax = T(2)
        mesh = UniformPeriodicMesh1D(xmin, xmax, N)
        cD_continuous = couple_continuosly(D, mesh)
        cD_central    = couple_discontinuosly(D, mesh)
        cD_plus       = couple_discontinuosly(D, mesh, Val(:plus))
        cD_minus      = couple_discontinuosly(D, mesh, Val(:minus))
        for cD in (cD_continuous, cD_central, cD_plus, cD_minus)
          print(devnull, cD)
          x = grid(cD)
          @test norm(cD * x.^0) < 100N * eps(T)
        end
        M = mass_matrix(cD_central)
        @test M ≈ mass_matrix(cD_plus)
        @test M ≈ mass_matrix(cD_minus)
        @test sum(M) ≈ xmax - xmin
        Dp = Matrix(cD_plus)
        Dm = Matrix(cD_minus)
        diss = M * (Dp - Dm)
        @test diss ≈ diss'
        @test maximum(eigvals(Symmetric(diss))) < 100N * eps(T)
        Dc = Matrix(cD_central)
        res = M * Dc + Dc' * M
        @test norm(res) < 100N * eps(T)
        res = M * Dp + Dm' * M
        @test norm(res) < 100N * eps(T)

        Mcont = mass_matrix(cD_continuous)
        @test sum(Mcont) ≈ xmax - xmin
        Dcont = Matrix(cD_continuous)
        res = Mcont * Dcont + Dcont' * Mcont
        @test norm(res) < 100N * eps(T)
      end
    end
  end
end


@testset "Coupled FD operators" begin
  for T in (Float32, Float64)
    for acc_order in 2:2:6
      for source in (MattssonNordström2004(),
                     MattssonAlmquistCarpenter2014Extended(),
                     MattssonAlmquistCarpenter2014Optimal())
        ymin = T(0)
        ymax = T(2)
        D = derivative_operator(source, 1, acc_order, ymin, ymax, 31)

        for N in 1:3
          xmin = T(-1)
          xmax = T(2)
          mesh = UniformMesh1D(xmin, xmax, N)
          cD_continuous = couple_continuosly(D, mesh)
          cD_central    = couple_discontinuosly(D, mesh)
          cD_plus       = couple_discontinuosly(D, mesh, Val(:plus))
          cD_minus      = couple_discontinuosly(D, mesh, Val(:minus))
          for cD in (cD_continuous, cD_central, cD_plus, cD_minus)
            print(devnull, cD)
            x = grid(cD)
            @test norm(cD * x.^0) < 100N * eps(T)

            u = sinpi.(x)
            v = copy(u)
            SummationByPartsOperators.scale_by_mass_matrix!(u, cD)
            SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, cD)
            @test u ≈ v
          end
          M = mass_matrix(cD_central)
          @test M ≈ mass_matrix(cD_plus)
          @test M ≈ mass_matrix(cD_minus)
          @test sum(M) ≈ xmax - xmin
          Dp = Matrix(cD_plus)
          Dm = Matrix(cD_minus)
          diss = M * (Dp - Dm)
          @test diss ≈ diss'
          @test maximum(eigvals(Symmetric(diss))) < 100N * eps(T)
          Dc = Matrix(cD_central)
          res = M * Dc + Dc' * M
          res[1, 1] += 1
          res[end, end] -= 1
          @test norm(res) < 100N * eps(T)
          res = M * Dp + Dm' * M
          res[1, 1] += 1
          res[end, end] -= 1
          @test norm(res) < 100N * eps(T)

          cD1 = couple_discontinuosly(cD_central, mesh)
          cD2 = couple_discontinuosly(D, UniformMesh1D(xmin, xmax, N^2))
          @test grid(cD1) ≈ grid(cD2)
          @test mass_matrix(cD1) ≈ mass_matrix(cD2)
          @test Matrix(cD1) ≈ Matrix(cD2)

          Mcont = mass_matrix(cD_continuous)
          @test sum(Mcont) ≈ xmax - xmin
          Dcont = Matrix(cD_continuous)
          res = Mcont * Dcont + Dcont' * Mcont
          res[1, 1] += 1
          res[end, end] -= 1
          @test norm(res) < 100N * eps(T)
        end

        for N in 1:3
          xmin = T(-1)
          xmax = T(2)
          mesh = UniformPeriodicMesh1D(xmin, xmax, N)
          cD_continuous = couple_continuosly(D, mesh)
          cD_central    = couple_discontinuosly(D, mesh)
          cD_plus       = couple_discontinuosly(D, mesh, Val(:plus))
          cD_minus      = couple_discontinuosly(D, mesh, Val(:minus))
          for cD in (cD_continuous, cD_central, cD_plus, cD_minus)
            print(devnull, cD)
            x = grid(cD)
            @test norm(cD * x.^0) < 100N * eps(T)
          end
          M = mass_matrix(cD_central)
          @test M ≈ mass_matrix(cD_plus)
          @test M ≈ mass_matrix(cD_minus)
          @test sum(M) ≈ xmax - xmin
          Dp = Matrix(cD_plus)
          Dm = Matrix(cD_minus)
          diss = M * (Dp - Dm)
          @test diss ≈ diss'
          @test maximum(eigvals(Symmetric(diss))) < 100N * eps(T)
          Dc = Matrix(cD_central)
          res = M * Dc + Dc' * M
          @test norm(res) < 100N * eps(T)
          res = M * Dp + Dm' * M
          @test norm(res) < 100N * eps(T)

          Mcont = mass_matrix(cD_continuous)
          @test sum(Mcont) ≈ xmax - xmin
          Dcont = Matrix(cD_continuous)
          res = Mcont * Dcont + Dcont' * Mcont
          @test norm(res) < 100N * eps(T)
        end
      end
    end
  end
end
