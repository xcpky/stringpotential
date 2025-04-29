using Base.Threads
using BenchmarkTools
using LinearAlgebra
using SpecialFunctions
using SphericalHarmonics
using Plots
using QuadGK
pyplot()
include("constants.jl")

const debug = false
const n_max = 30
const r_1 = 0.02 * 5.068  # GeV^-1
const r_n_max = 30 * 5.068  # GeV^-1
const c_T = (-1) * 0.19732^2 / (1 / 2 * 4.18) * 5.068

function quadgauss(f, x::T, w::T) where {T<:Vector{Float64}}
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    @inbounds @simd for i in eachindex(x)
        res += f(x[i]) * w[i]
    end
    return res
end



# Double factorial function
function doublefactorial(n::Int)
    n == 0 || n == 1 ? 1 : n * doublefactorial(n - 2)
end

# Radial coordinate functions
r_n(n) = r_1 * (r_n_max / r_1)^((n - 1) / (n_max - 1))
nu_n(n) = 1 / r_n(n)^2

# Normalization functions
N_nl(n, l) = sqrt(2^(l + 2) * (2 * nu_n(n))^(l + 3 / 2) / (sqrt(π) * doublefactorial(2l + 1)))

N_n_n_prime(n, n_prime, l) = (2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime)))^(l + 3 / 2)

T_n_n_prime(n, n_prime, l) = c_T * (2l + 3) * nu_n(n) * nu_n(n_prime) / (nu_n(n) + nu_n(n_prime)) * N_n_n_prime(n, n_prime, l)

V_r_n_n_prime(n, n_prime, l) = 1 / sqrt(2 * (nu_n(n) + nu_n(n_prime))) *
                               doublefactorial(2l + 2) / doublefactorial(2l + 1) * N_n_n_prime(n, n_prime, l)

V_1_r_n_n_prime(n, n_prime, l) = 2 / √π * 2^l * factorial(l) / doublefactorial(2l + 1) *
                                 sqrt(nu_n(n) + nu_n(n_prime)) * N_n_n_prime(n, n_prime, l)

# Modified normalization functions
N_nl_tilde(n, l) = sqrt(2^(l + 3) * (2 * nu_n(n))^(l + 5 / 2) / (sqrt(π) * doublefactorial(2l + 3)))

N_n_n_prime_tilde(n, n_prime, l) = (2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime)))^(l + 5 / 2)

V_r_n_n_prime_tilde(n, n_prime, l) = 1 / sqrt(2 * (nu_n(n) + nu_n(n_prime))) *
                                     doublefactorial(2l + 4) / doublefactorial(2l + 3) * N_n_n_prime_tilde(n, n_prime, l)

# Potential functions
V_e_mu_n_n_prime(n, n_prime, l, μ) = (2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime) + μ))^(l + 3 / 2)

function V_e_mu_lin_n_n_prime_l_0(n, n_prime, μ)
    term1 = 1 / (2 * (nu_n(n) + nu_n(n_prime)))
    term2 = sqrt(π) * μ * erfcx(μ / (2 * sqrt(nu_n(n) + nu_n(n_prime)))) / (4 * (nu_n(n) + nu_n(n_prime))^(3 / 2))
    sqrt(4 * (2 * nu_n(n))^(3 / 2) / sqrt(π)) * sqrt(4 * (2 * nu_n(n_prime))^(3 / 2) / sqrt(π)) * (term1 - term2)
end

V_deuteron_n_n_prime(n, n_prime) = -626.885 * V_e_mu_lin_n_n_prime_l_0(n, n_prime, 1.55) +
                                   1438.72 * V_e_mu_lin_n_n_prime_l_0(n, n_prime, 3.11)

# Set up and solve eigenvalue problem
l = 0
Lam = 0.0
scale_factor = 1.0
V_0_fit = 0.4

V_0_l = -0.08858455651140933
sigma_l = 0.0199179973550142
alpha_l = 0.476814032326273
r_0_l = 15
V_0 = V_0_l / a_l
sigma = sigma_l / (a_l^2)
alpha = alpha_l
r_0 = r_0_l * a_l
# Initialize matrices
H_Cornell = zeros(n_max, n_max)
N_n_n = zeros(n_max, n_max)

for i in 1:n_max
    n = i
    for j in 1:n_max
        n_prime = j
        H_ij = -T_n_n_prime(n, n_prime, l) + V_0_fit * N_n_n_prime(n, n_prime, l) +
               sigma * V_r_n_n_prime(n, n_prime, l) - alpha * V_1_r_n_n_prime(n, n_prime, l)
        N_nn = N_n_n_prime(n, n_prime, l)

        H_Cornell[i, j] = H_ij
        N_n_n[i, j] = N_nn
    end
end

# Solve generalized eigenvalue problem
E_solution = eigen(H_Cornell, N_n_n).values
c_solution = eigen(H_Cornell, N_n_n).vectors
for i in 1:n_max
    c_solution[:, i] ./= sqrt(c_solution[:, i]' * N_n_n * c_solution[:, i])
end

# Wavefunction definitions
function psi_nlm(r, n, l, m, θ, φ)
    ψ = 0.0
    for i in 1:n_max
        n_prime = i
        ψ += c_solution[i, n] * N_nl(n_prime, l) * r^l * exp(-nu_n(n_prime) * r^2) *
             SphericalHarmonics.sphericalharmonic(θ, φ, l=l, m=m)
    end
    ψ
end

function psi_nlm_tilde(r, n, l, m, θ, φ)
    ψ = 0.0
    for i in 1:n_max
        n_prime = i
        ψ += c_solution[i, n] * N_nl_tilde(n_prime, l) * r^(l + 1) *
             exp(-nu_n(n_prime) * r^2) * SphericalHarmonics.sphericalharmonic(θ, φ, l=l, m=m)
    end
    ψ
end

# Momentum space wavefunction
N_nl_p(n_prime) = sqrt(4 * sqrt(2) / (9 * π^(7 / 2)) * (nu_n(n_prime))^(5 / 2))

function psi_nlm_cs(r, n, l, m, θ, ϕ)
    if psi_nlm(0.0001, n, l, m, θ, ϕ) < 0
        return -psi_nlm(r, n, l, m, θ, ϕ)
    else
        return psi_nlm(r, n, l, m, θ, ϕ)
    end
end

function psi_nlm_p(p, n)
    ψ_p = 0.0 + 0.0im
    for i in 1:n_max
        n_prime = i
        ψ_p += c_solution[i, n] * N_nl_p(n_prime) * sqrt(3) * π * im / (4 * nu_n(n_prime)^(5 / 2)) *
               p * exp(-p^2 / (4 * nu_n(n_prime)))
    end
    real(psi_nlm(1e-4, n, 1, 0, 0, 0)) < 0.0 ? -ψ_p : ψ_p
end

if debug

    r = 0:0.0001:7
    p = plot(r, r .* real.(psi_nlm_cs.(r, 1, 0, 0, 0, 0)) ./ real(SphericalHarmonics.sphericalharmonic(0, 0, 0, 0)), grid=true, dpi=300)
    for i in 2:7
        plot!(p, r, r .* real.(psi_nlm_cs.(r, i, 0, 0, 0, 0)) ./ real(SphericalHarmonics.sphericalharmonic(0, 0, 0, 0)))
    end
    savefig(p, "plt.png")
end

function PIgauss(x, a1, a2, mu1, mu2, sigma1, sigma2)
    a1 * exp(-(x - mu1)^2 / (2 * sigma1^2)) + a2 * exp(-(x - mu2)^2 / (2 * sigma2^2))
end

function psi_3P1_momentum(p, n, Ngauss, Λ)
    x, w = gauss(Ngauss, 0, Λ)
    integrand1(x) = exp(1im * x * p) * x / p * psi_nlm_cs(x, n, 1, 0, 0, 0)
    integrand2(x) = -exp(-1im * x * p) * x / p * psi_nlm_cs(x, n, 1, 0, 0, 0)
    result1 = -im * 2 * π * (quadgauss(integrand1, x, w) + quadgauss(integrand2, x, w))
    result2 = real(result1)^2 + imag(result1)^2
    return result1, result2
end

function psi_S_1_3_momentum(p, n, Ngauss, Λ)
    x, w = gauss(Ngauss, 0, Λ)
    integrand1(x) = exp(im * x * p) * x / p * psi_nlm_cs(x, n, 0, 0, 0, 0)
    integrand2(x) = -exp(-im * x * p) * x / p * psi_nlm_cs(x, n, 0, 0, 0, 0)
    result1 = -im * 2 * π * (quadgauss(integrand1, x, w) + quadgauss(integrand2, x, w))
    result2 = real(result1)^2 + imag(result1)^2
    return result1, result2
end

E_sample = range(delta[1]-1, delta[2] + 0.4, 1000)
Nq = 40
Ntower = 30
xq, wq = gauss(Nq)
psi_mat_q1 = zeros((Ntower * length(E_sample), Nq + 1))
psi_mat_q2 = zeros((Ntower * length(E_sample), Nq + 1))
@inline function xsqrt(x)
    imag(x) >= 0 ? sqrt(x + 0im) : -sqrt(x - 0im)
end
@inline function q0(E)
    dE = -(delta .- E)
    x0 = xsqrt.(2 * mu .* dE)
    return x0
end

@time for i in 1:Ntower
    for k in 1:psi_mat_q1.size[2]-1
        psi_mat_q1[i, k] = real(psi_3P1_momentum(xq[k], i, Nq, 20)[1])
        psi_mat_q2[i, k] = psi_mat_q1[i, k]
    end
    @threads for l in eachindex(E_sample)
        psi_mat_q1[i+Ntower*(l-1), Nq+1] = real(psi_3P1_momentum(q0(E_sample[l])[1], i, Nq, 20)[1])
        psi_mat_q2[i+Ntower*(l-1), Nq+1] = real(psi_3P1_momentum(q0(E_sample[l])[2], i, Nq, 20)[1])
        psi_mat_q1[i+Ntower*(l-1), 1:end-1] = psi_mat_q1[i, 1:end-1]
        psi_mat_q2[i+Ntower*(l-1), 1:end-1] = psi_mat_q2[i, 1:end-1]
    end
end

# for i in 1:Ntower
#     for k in 1:psi_mat_q1.size[2]-1
#         psi_mat_q1[i, k] = real(psi_nlm_p(xq[k], i))
#         psi_mat_q2[i, k] = psi_mat_q1[i, k]
#     end
#     for l in eachindex(E_sample)
#         psi_mat_q1[i+Ntower*(l-1), Nq+1] = imag(psi_nlm_p(q0(E_sample[l])[1], i))
#         psi_mat_q2[i+Ntower*(l-1), Nq+1] = imag(psi_nlm_p(q0(E_sample[l])[2], i))
#         psi_mat_q1[i+Ntower*(l-1), 1:end-1] = psi_mat_q1[i, 1:end-1]
#         psi_mat_q2[i+Ntower*(l-1), 1:end-1] = psi_mat_q2[i, 1:end-1]
#     end
# end
