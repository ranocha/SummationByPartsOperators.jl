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


# 6th order dissipation
D3 = zeros(Int, N, N)
for i in 1:N
    D3[i,i] = -3
end
for i in 1:N-1
    D3[i+1,i] = 3
    D3[i,i+1] = 1
end
for i in 1:N-2
    D3[i+2,i] = -1
end
boundary_width = 2
for i in 1:boundary_width
    D3[i,:] = D3[boundary_width+1,:]
end
for i in 1:boundary_width-1
    D3[end-i+1,:] = D3[end-boundary_width+1,:]
end
display(D3)

B3 = eye(Int, N); B3[1,1] = B3[2,2]= B3[end,end] = 0
display(B3)

display(D3' * B * D3)
display(D3' * B3 * D3)


# 8th order dissipation
D4 = zeros(Int, N, N)
for i in 1:N
    D4[i,i] = 6
end
for i in 1:N-1
    D4[i+1,i] = -4
    D4[i,i+1] = -4
end
for i in 1:N-2
    D4[i+2,i] = 1
    D4[i,i+2] = 1
end
boundary_width = 2
for i in 1:boundary_width
    D4[i,:] = D4[boundary_width+1,:]
    D4[end-i+1,:] = D4[end-boundary_width,:]
end
display(D4)

B4 = eye(Int, N); B4[1,1] = B4[2,2]= B4[end,end] = B4[end-1,end-1] = 0
display(B4)

display(D4' * B * D4)
display(D4' * B4 * D4)
