using DynamicPolynomials
using NCTSSOS

n = 3
@ncpolyvar x[1:3]
f = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3.0x[2]^2 - 2x[1]*x[2]*x[1] + 2x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] +
6x[3]^2 + 9x[2]^2*x[3] + 9x[3]*x[2]^2 - 54x[3]*x[2]*x[3] + 142x[3]*x[2]^2*x[3]
pop = [f]
opt,data = ncpop(pop, x, 2, CS=false, TS="block", newton=true, Gram=true)
# opt = 0 

for i = 1:n
    push!(pop, 1 - x[i]^2)
    push!(pop, x[i] - 1/3)
end
opt,data = ncpop(pop, x, 2, TS="block", Gram=true)
# opt = 0.9975308
opt,data = ncpop(data, TS="block", Gram=true)

n = 2
@ncpolyvar x[1:2]
f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2 + x[1]*x[2]*x[1]*x[2] + x[2]*x[1]*x[2]*x[1] +
x[1]^3*x[2] + x[2]*x[1]^3 + x[1]*x[2]^3 + x[2]^3*x[1]
g1 = 1 - x[1]^2
g2 = 1 - x[2]^2
pop = [f, g1, g2]
opt,data = ncpop(pop, x, 2, TS="block")
# opt = -2.051110

f = x[2]*x[1] + x[2]*x[1]^2*x[2] - 2x[2]*x[1]^2 + x[2]*x[1]*x[2]*x[1] - 2x[2]*x[1]*x[2] +
x[1]*x[2] + x[1]*x[2]*x[1]*x[2] - 2x[1]^2*x[2] + x[2]*x[1]^2*x[2] - 2x[2]*x[1]*x[2]
pop = [f, 4-x[1]^2, 4-x[2]^2]
opt,data = ncpop(pop, x, 2, TS="block", obj="trace")
# opt = -8

@ncpolyvar x[1:6]
f = x[1]*(x[4] + x[5] + x[6]) + x[2]*(x[4] + x[5] - x[6]) + x[3]*(x[4] - x[5]) - x[1] - 2x[4] - x[5]
opt,data = ncpop([-f], x, 2, partition=3, constraint="projection", TS="block", QUIET=true, obj="eigen")
# opt = -0.2590717
opt,data = ncpop([-f], x, 3, soc=true, partition=3, constraint="projection", CS=false, TS="block", QUIET=true, obj="eigen")
# opt = -0.25087557

# Broyden banded polynomial
n = 10
@ncpolyvar x[1:n]
f = 0.0
for i = 1:n
    jset = max(1, i-5) : min(n, i+1)
    jset = setdiff(jset, i)
    g = sum(x[j] + x[j]^2 for j in jset)
    f += (2*x[i] + 5*x[i]^3 + 1 - g)^2
end
pop = [f]
for i = 1:n
    push!(pop, 1 - x[i]^2)
    push!(pop, x[i] - 1/3)
end
opt,data = ncpop(pop, x, 3, TS="MD")
# opt = 3.011288


# Chained singular polynomial
n = 10
f = 0.0
@ncpolyvar x[1:n]
for i = 1:2:n-3
    f += (x[i] + 10*x[i+1])^2 + 5*(x[i+2] - x[i+3])^2 + (x[i+1] - 2*x[i+2])^4 + 10*(x[i] - 10*x[i+3])^4
end
pop = [f]
for i = 1:n
    push!(pop, 1 - x[i]^2)
    push!(pop, x[i] - 1/3)
end
opt,data = ncpop(pop, x, 2, TS="MD")
# opt = 81.18341
