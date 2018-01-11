# order 2
using SymPy

N = 11
b = [symbols("b_$i", real=true) for i in 1:N]

# order 2
h = [1//2, 1, 1, 1]
H = Diagonal(h)

s = [-3//2, 2, -1//2, 0]

m1 = [(b[1]+b[2])/2, -(b[1]+b[2])/2, 0, 0]
m2 = [-(b[1]+b[2])/2, b[1]/2+b[2]+b[3]/2, -(b[2]+b[3])/2, 0]
m3 = [0, -(b[2]+b[3])/2, b[2]/2+b[3]+b[4]/2, ]
M = vcat(m1, m2, m3)'

D2 = -H \ M; D2[1,:] -= b[1]*s
