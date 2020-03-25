using LinearAlgebra

# order 2
Qp = [
  -1//4 5//4 -1//2 0 0 0
  -1//4 -5//4 2 -1//2 0 0
  0 0 -3//2 2 -1//2 0
  0 0 0 -3//2 2 -1//2
  0 0 0 0 -5//4 5//4
  0 0 0 0 -1//4 -1//4
]
Qm = Matrix(-Qp')
B = zero(Qp); B[1,1] = -1; B[end,end] = 1
M = Diagonal([1//4, 5//4, 1, 1, 5//4, 1//4])
Dp = M \ (Qp + B / 2)
Dm = M \ (Qm + B / 2)

# order 3
