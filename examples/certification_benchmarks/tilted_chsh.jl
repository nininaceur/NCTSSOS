import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using CertifiedQuantumBounds
using DynamicPolynomials
using NCTSSOS
import Polynomials: roots
import Polynomials
const PPoly = Polynomials.Polynomial

@ncpolyvar X[1:2]
@ncpolyvar Y[1:2]

CHSHt(alpha, beta) = X[1]*Y[1] +  X[1]*Y[2]  + X[2]*Y[1]  - X[2]*Y[2] + alpha*X[1] + beta*Y[1]

function max_violation(alpha, beta)
    tau0 = alpha^6*(27*beta^2 - 8) + alpha^4*(-54*beta^4 + 48*beta^2 + 32.0) + alpha^2*(27*beta^6 + 48*beta^4 - 400*beta^2 + 128) - 8*beta^6 + 32*beta^4 + 128*beta^2 - 512.0
    tau1 = -168*alpha^3*beta^3 + 60*alpha^5*beta + 160*alpha^3*beta + 60*alpha*beta^5 + 160*alpha*beta^3 - 576*alpha*beta
    tau2 = -6*alpha^4*beta^2 - 6*alpha^2*beta^4 - 64*alpha^2*beta^2 + 20*alpha^4 + 96*alpha^2 + 20*beta^4 + 96*beta^2 + 320
    tau3 =  8*alpha^3*beta^3 - 24*alpha^3*beta - 24*alpha*beta^3 + 96*alpha*beta
    tau4 = 11*alpha^2*beta^2 - 16*alpha^2 - 16*beta^2 - 64
    tau5 = 4*alpha*beta
    tau6 = 4.0

    coeffs = [tau0, tau1, tau2, tau3, tau4, tau5, tau6]    # ascending-order
    p = DynamicPolynomials.Polynomial(coeffs)
    lambda      = roots(p)
    reallambda  = real.(lambda[abs.(imag.(lambda)) .< 1e-10])
    return maximum(reallambda)
end

function max_violation_sym(alpha::Real)
    alpha = float(alpha)

    # coefficients of lambda⁰…lambda⁴
    c0 =  5*alpha^6 - 21*alpha^4 + 16*alpha^2 - 32
    c1 =  2*alpha^6 -    alpha^4 - 20*alpha^2 - 32
    c2 =  11/4*alpha^4 - 12*alpha^2 - 4
    c3 =  4      -    alpha^2
    c4 =  1.0

    p = PPoly([c0, c1, c2, c3, c4])
    lambdas = roots(p)
    reallambda = real.(lambdas[abs.(imag.(lambdas)) .< 1e-8])
    return maximum(reallambda)
end

function bounds_table(alpha, r_max::Int)
    tbl = Array{NTuple{8,Float64}}(undef, 1, r_max)   # (new, old, exact)
    for r in 1:r_max
        println("\nRelaxation order: ", r)
        newB, oldB, Bdiff, eigproj, eigraw, diffraw, glength = rational_certificate(
                        -CHSHt(alpha,alpha), [], [], [X;Y], r;
                        partition = 2,
                        constraint = "unipotent",
                        QUIET = false,
                        QUIETTS = true,
                        tol = 1e-25)
        tbl[r] = (newB, oldB, max_violation_sym(alpha), Bdiff, eigproj, eigraw, diffraw, glength)
    end
    return tbl
end

alph = 0.999
r_max = 10

tbl = bounds_table(alph, r_max)
