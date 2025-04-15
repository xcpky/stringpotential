module SchrodingerEquation
using QuadGK
using LinearAlgebra
using GSL

# Export public symbols
export solve, Solver, M_0, M_l

# Define a type
"""
    Solver(potential::Function, boundary_conditions::Vector{Float64}, grid_spacing::Float64=0.01)

A solver for the Schrödinger equation.

# Arguments
- `potential::Function`: The potential energy function `V(x)`.
- `boundary_conditions::Vector{Float64}`: Boundary conditions for the problem.
- `grid_spacing::Float64=0.01`: Spacing for the spatial grid (default: 0.01).
"""
struct Solver
    Ngauss::Int
    p::Vector{Float64}
    w::Vector{Float64}
    y::Matrix{Float64}
    KinCoeffi::Float64
    LinearCoeffi::Float64
    ReciproCoeffi::Float64
    S::Array{Array{Float64}}
    dif1::Matrix{Float64}
    dif2::Matrix{Float64}
    LegendreMat::Array{Float64}
    LegendreDerivMat::Array{Float64}
    Wderiv::Array{Float64}
    Qmat::Array{Float64}
    Q0::Array{Float64}

    # Constructor
    function Solver(Scaling::Float64, KinCoeffi::Float64, LinearCoeffi::Float64, ReciproCoeffi::Float64, Ngauss::Int64=800, Nlag::Int64=5)
        # Initialize energy levels (example)
        p, w = gauss(Ngauss)
        @. p = Scaling * ((1 + p) / (1 - p))
        @. w = Scaling * (2 * w / (1 - p)^2)
        y = @. (p^2 + p'^2) / 2 / (p * p')
        s1 = zeros(Float64, Ngauss)
        s2 = zeros(Float64, Ngauss)
        s3 = zeros(Float64, Ngauss)
        @inbounds for i in 1:Ngauss
            @inbounds for j in 1:Ngauss
                if i == j
                    continue
                end
                s1[i] += w[j] * 2 * p[j]^2 / (p[j]^2 - p[i]^2)^2
                s2[i] += w[j] / (p[j]^2 - p[i]^2)
                yij = (p[i]^2 + p[j]^2) / 2 / p[i] / p[j]
                s3[i] += w[j] / p[j] / 2 * log((yij + 1) / (yij - 1))
            end
        end
        dif1, dif2 = Dmat(Ngauss, Nlag, p)
        LegendreMat, LegendreDerivMat, Qmat, Q0, Wderiv = legendre(Ngauss, y)
        new(Ngauss, p, w, y, KinCoeffi, LinearCoeffi, ReciproCoeffi, [s1, s2, s3], dif1, dif2, LegendreMat, LegendreDerivMat, Wderiv, Qmat, Q0)
    end
end

function solve(solver::Solver, l::Int)
    M = if l == 0
        M_0(solver)
    else
        M_l(solver, l)
    end
    eig_result = eigen(Hermitian(M))  # Use Hermitian wrapper if applicable

    # Extract and sort eigenvalues (ascending order)
    sorted_indices = sortperm(eig_result.values)
    eigenvalues = eig_result.values[sorted_indices]
    eigenvectors = eig_result.vectors[:, sorted_indices]
    return eigenvalues, eigenvectors
end

function W1(l)
    sum([1 / m for m in 1:l])
end

function M_0(solver::Solver)
    Ngauss = solver.Ngauss
    KinCoeffi = solver.KinCoeffi
    LinearCoeffi = solver.LinearCoeffi
    ReciproCoeffi = solver.ReciproCoeffi
    p = solver.p
    S = solver.S
    w = solver.w
    Wderiv = solver.Wderiv
    dif1 = solver.dif1
    dif2 = solver.dif2
    ql = solver.Qmat
    q0 = solver.Q0
    M = Matrix{Float64}(undef, (Ngauss, Ngauss))
    Pl = solver.LegendreMat
    Plderiv = solver.LegendreDerivMat
    for i in 1:Ngauss
        for j in 1:Ngauss
            if i == j
                M[i, j] = KinCoeffi * p[i]^2 + ReciproCoeffi * p[i] / pi * (S[3][i] - pi^2 / 2) + 2 * LinearCoeffi / pi * S[1][i]
            else
                M[i, j] = -ReciproCoeffi / pi * w[j] * p[j] / p[i] * q0[i, j]
                M[i, j] -= 2 * LinearCoeffi / pi * w[j] * (2 * p[j]^2 / (p[j]^2 - p[i]^2)^2)
            end

            M[i, j] += 2 * LinearCoeffi / pi * ((p[i] * S[2][i] - 3 * w[i] / 4 / p[i]) * dif1[i, j] - w[i] / 4 * dif2[i, j])
        end
    end
    return M
end

function M_l(solver::Solver, l::Int)
    Ngauss = solver.Ngauss
    KinCoeffi = solver.KinCoeffi
    LinearCoeffi = solver.LinearCoeffi
    ReciproCoeffi = solver.ReciproCoeffi
    p = solver.p
    S = solver.S
    w = solver.w
    Wderiv = solver.Wderiv
    dif1 = solver.dif1
    dif2 = solver.dif2
    ql = solver.Qmat
    q0 = solver.Q0
    M = Matrix{Float64}(undef, (Ngauss, Ngauss))
    Pl = solver.LegendreMat
    Plderiv = solver.LegendreDerivMat
    for i in 1:Ngauss
        for j in 1:Ngauss
            if i == j
                M[i, j] = KinCoeffi * p[i]^2 + ReciproCoeffi * p[i] / pi * (S[3][i] - pi^2 / 2 + w[i] * W1(l) / p[i]) + 2 * LinearCoeffi / pi * (S[1][i] + Plderiv[l, i, i] / 2 / p[i] * (pi^2 / 2 - S[3][i] - w[i] / 2 / p[i]))
            end
            if i != j
                M[i, j] = -ReciproCoeffi / pi * w[j] * p[j] / p[i] * ql[l, i, j]
                M[i, j] -= 2 * LinearCoeffi / pi * w[j] * (2 * p[j]^2 / (p[j]^2 - p[i]^2)^2 * Pl[l, i, j] - q0[i, j] / 2 / p[i]^2 * Plderiv[l, i, j])
            end
                M[i, j] += 2 * LinearCoeffi / pi * ((p[i] * S[2][i] - 3 * w[i] / 4 / p[i]) * dif1[i, j] - w[i] / 4 * dif2[i, j] - w[j] * Wderiv[l, i, j] / 2 / p[i]^2)
        end
    end
    return M
end


function Dmat(Ngauss::Int64, Nlag::Int64, p::Vector{Float64})
    dif1 = Matrix{Float64}(undef, Ngauss, Ngauss)
    dif2 = Matrix{Float64}(undef, Ngauss, Ngauss)
    for j in 1:Ngauss
        for i in 1:Ngauss
            start = min(Ngauss - Nlag + 1, max(1, i - Int(floor(Nlag // 2))))
            range = start:start+Nlag-1
            mask = setdiff(range, [i, j])
            if i != j
                local temp = (p[i] .- p[mask]) ./ (p[j] .- p[mask])
                dif1[i, j] = 1 / (p[j] - p[i]) * prod(temp)
                dif2[i, j] = 0
                for m in mask
                    mas = setdiff(mask, m)
                    local temp = (p[i] .- p[mas]) ./ (p[j] .- p[mas])
                    dif2[i, j] += 1 / (p[j] - p[m]) * prod(temp)
                end
                # for l in mask
                #     mas = setdiff(mask, l)
                #     dif2[i, j] += 1 / (p[j] - p[l]) * prod((p[mas] .- p[i]) ./ (p[mas] .- p[j]))
                # end
                dif2[i, j] *= 2 / (p[j] - p[i])
            else
                local temp = 1 ./ (p[j] .- p[mask])
                dif1[i, i] = sum(temp)
                dif2[i, i] = 0
                for m in mask
                    mas = setdiff(mask, m)
                    local temp = 1 ./ (p[j] .- p[mas])
                    dif2[i, i] += 1 / (p[j] - p[m]) * sum(temp)
                end
            end
        end
    end
    return dif1, dif2
end

function legendre(Ngauss, y)
    # Initialize arrays with the 3-length dimension first
    legendre_result = Array{Float64}(undef, (3, Ngauss, Ngauss))
    legendre_deriv_result = Array{Float64}(undef, (3, Ngauss, Ngauss))
    q = Array{Float64}(undef, (3, Ngauss, Ngauss))
    q_0 = Matrix{Float64}(undef, (Ngauss, Ngauss))
    w = Array{Float64}(undef, (3, Ngauss, Ngauss))

    # Populate the arrays
    legendre_result[1, :, :] = y
    legendre_result[2, :, :] = GSL.sf_legendre_P2.(y)
    legendre_result[3, :, :] = GSL.sf_legendre_P3.(y)

    legendre_deriv_result[1, :, :] .= 1
    legendre_deriv_result[2, :, :] = sf_legendre_deriv_P2.(y)
    legendre_deriv_result[3, :, :] = sf_legendre_deriv_P3.(y)

    for i in 1:Ngauss
        for j in 1:i-1
            q[1, i, j] = GSL.sf_legendre_Q1(y[i, j])
            q[2, i, j] = GSL.sf_legendre_Ql(2, y[i, j])
            q[3, i, j] = GSL.sf_legendre_Ql(3, y[i, j])
            q_0[i, j] = GSL.sf_legendre_Q0(y[i, j])
        end
        q[:, i, i] .= 0
        q_0[i, i] = 0
    end

    w[1, :, :] = zeros((Ngauss, Ngauss))
    for l in 2:3
        for i in 1:Ngauss, j in 1:i
            w[l, i, j] = legendre_deriv_result[l-1, i, j]
            for m in 2:l-1
                w[l, i, j] += (legendre_deriv_result[l-m, i, j] * legendre_result[m-1, i, j] + legendre_result[l-m, i, j] * legendre_deriv_result[m-1, i, j]) / m
            end
            if l > 1
                w[l, i, j] += legendre_deriv_result[l-1, i, j] / l
            end
        end
    end

    # Ensure symmetry for each sub-array (i in 1:3)
    for i in 1:3
        legendre_result[i, :, :] = Symmetric(legendre_result[i, :, :], :L)
        legendre_deriv_result[i, :, :] = Symmetric(legendre_deriv_result[i, :, :], :L)
        q[i, :, :] = Symmetric(q[i, :, :], :L)
        w[i, :, :] = Symmetric(w[i, :, :], :L)
    end

    return legendre_result, legendre_deriv_result, q, Symmetric(q_0, :L), w
end

@inline function sf_legendre_deriv_P2(x::Float64)
    3x
end

@inline function sf_legendre_deriv_P3(x::Float64)
    (15x^2 - 3) / 2
end

@generated function dfactorial(::Val{n}) where {n}
    # Compute the result at compile time
    result = 1
    current = n
    while current > 0
        result *= current
        current -= 2
    end
    return :($result)
end

# Helper for cleaner syntax
dfactorial(n::Int) = dfactorial(Val{n}())

function quadgauss(f, x::Vector{Float64}, w::Vector{Float64})
    # Type-stable initialization
    T = typeof(f(x[1]) * w[1])
    res = zero(T)

    # SIMD-accelerated, bounds-check-free loop
    @inbounds @simd for i in eachindex(x)
        res += f(x[i]) * w[i]
    end

    return res
end

function Q_0(x)
    1 / 2 * log((x + 1) / (x - 1))
end

# Include additional files (optional)
end
