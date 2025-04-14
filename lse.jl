# include("hadronmasses.jl")
# include("constants.jl")
include("stringpotential.jl")
using LinearAlgebra
using QuadGK

function quadgauss(f, x::T, w::T) where {T<:Vector{Float64}}
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    for i in eachindex(x)
        res += f(x[i]) * w[i]
    end
    return res
end

function gmat(Λ, E, Ngauss)
    nodes, weights = gauss(Ngauss, 0, Λ)
    dim = 2 * Ngauss + 2
    mat = zeros(ComplexF64, dim, dim)
    dE = -(delta .- E)
    x0 = xsqrt.(2 * mu .* dE)
    if E >= delta[1]
        g1(x) = 1 / (dE[1] - x^2 / 2 / mu[1] + im * ϵ)
        integ = quadgauss(g1, nodes, weights)
        mat[Ngauss+1, Ngauss+1] = 1 / (2 * π^2) * x0[1] * (mu[1] * (-im * π + log((Λ + x0[1]) / (Λ - x0[1]))) - x0[1] * integ)
    else
        mat[dim, dim] = 0 + im * 0
    end
    if E >= delta[2]
        g2(x) = 1 / (dE[2] - x^2 / 2 / mu[2] + im * ϵ)
        integ = quadgauss(g2, nodes, weights)
        mat[dim, dim] = 1 / (2 * π^2) * x0[2] * (mu[2] * (-im * π + log((Λ + x0[2]) / (Λ - x0[2]))) - x0[2] * integ)
    else
        mat[dim, dim] = 0 + im * 0
    end
    for i = 1:Ngauss
        mat[i, i] = 1 / (2 * π^2) * nodes[i]^2 * weights[i] / (dE[1] - nodes[i]^2 / 2 / mu[1] + im * ϵ)
        idx = i + Ngauss + 1
        mat[idx, idx] = 1 / (2 * π^2) * nodes[i]^2 * weights[i] / (dE[2] - nodes[i]^2 / 2 / mu[2] + im * ϵ)
    end
    return mat
end

function vmat(Λ, E, Ngauss, g_C)
    nodes, _ = gauss(Ngauss, 0, Λ)
    p::Array{Vector{ComplexF64}} = [copy(nodes), copy(nodes)]
    dim = 2 * Ngauss + 2
    mat = zeros(ComplexF64, dim, dim)
    dE = -(delta .- E)
    x0 = xsqrt.(2 * mu .* dE)
    append!(p[1], x0[1])
    append!(p[2], x0[2])
    for i = 1:Ngauss+1
        for j = 1:Ngauss+1
            mat[i, j] = V_OME_11(E, p[1][i], p[1][j]) + Ctct_11(g_C)
            mat[i+Ngauss+1, j] = V_OME_21(E, p[2][i], p[1][j]) + Ctct_12(g_C)
            mat[i, j+Ngauss+1] = V_OME_12(E, p[1][i], p[2][j]) + Ctct_12(g_C)
            mat[i+Ngauss+1, j+Ngauss+1] = V_OME_22(E, p[2][i], p[2][j]) + Ctct_22(g_C)
        end
    end
    return mat
end

function tmat(Λ, E, Ngauss)
    G = gmat(Λ, E, Ngauss)
    V = vmat(Λ, E, Ngauss, 1)
    return inv(I - V * G) * V
end
