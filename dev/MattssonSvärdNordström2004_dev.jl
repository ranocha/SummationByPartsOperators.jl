using SymPy

N = 15

b = [symbols("b_$i", real=true) for i in 1:N]
B = diagm(b)


# 2nd order dissipation
D1 = zeros(Int, N, N)
for i in 1:N
    D1[i,i] = 1
end
for i in 1:N-1
    D1[i+1,i] = -1
end
D1[1,2] = 1
display(D1)

B1 = eye(Int, N); B1[1,1] = 0
display(B1)

display(D1' * B * D1)
display(D1' * B1 * D1)
