# order 2
h = [5//12, 7//6, 11//12, 1, 1]
e = [1//1, 0, 0, 0, 0]
H = Diagonal(h)

q1 = [0 7//12 -1//12 0 0]
q2 = [-q1[2] 0 7//12 0 0]
q3 = [-q1[3] -q2[3] 0 1//2 0]
q4 = [0 0 -1//2 0 1//2]
q5 = [0 0 0 -1//2 0]
Q = vcat(q1, q2, q3, q4, q5)

D1 = H \ (Q - 1//2*e*e'); display(D1)

x=1//1:length(h)
ks = 0:5
for k in ks
    u=x.^k
    println("k = ", k)
    println(D1*u - k*x.^(k-1))
end
