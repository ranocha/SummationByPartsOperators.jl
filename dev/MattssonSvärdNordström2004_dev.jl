using SymPy

N = 11

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


# 4th order dissipation
D2 = zeros(Int, N, N)
for i in 1:N
    D2[i,i] = -2
end
for i in 1:N-1
    D2[i+1,i] = 1
    D2[i,i+1] = 1
end
boundary_width = 1
for i in 1:boundary_width
    D2[i,:] = D2[boundary_width+1,:]
    D2[end-i+1,:] = D2[end-boundary_width,:]
end
display(D2)

B2 = eye(Int, N); B2[1,1] = B2[end,end] = 0
display(B2)

display(D2' * B * D2)
display(D2' * B2 * D2)
