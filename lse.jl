# include("hadronmasses.jl")
# include("constants.jl")
include("stringpotential.jl")
using LinearAlgebra
using QuadGK

function quadgauss(f, x::T, w::T) where {T<:Vector{Float64}}
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    @inbounds @simd for i in eachindex(x)
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
    E = E + m11 + m12
    for i = 1:Ngauss+1
        for j = 1:Ngauss+1
            # if i == 1 && j == 1
            #     println("V_OME_11($(E), $(p[1][i]), $(p[1][j])) + Ctct_11($(g_C))")
            #     println(V_OME_11(E, p[1][i], p[1][j]) + Ctct_11(g_C))
            #     println("V_OME_22($(E), $(p[2][i]), $(p[2][j])) + Ctct_22($(g_C))")
            #     println(V_OME_22(E, p[2][i], p[2][j]) + Ctct_22(g_C))
            # end
            # mat[i, j] = V_OME_11(E, p[1][i], p[1][j]) #+ Ctct_11(g_C)
            # mat[i+Ngauss+1, j] = V_OME_21(E, p[2][i], p[1][j]) #+ Ctct_12(g_C)
            # mat[i, j+Ngauss+1] = V_OME_12(E, p[1][i], p[2][j]) #+ Ctct_12(g_C)
            # mat[i+Ngauss+1, j+Ngauss+1] = V_OME_22(E, p[2][i], p[2][j]) #+ Ctct_22(g_C)
            mat[i, j] = V11(E, p[1][i], p[1][j]) #+ Ctct_11(g_C)
            mat[i+Ngauss+1, j] = V21(E, p[2][i], p[1][j]) #+ Ctct_12(g_C)
            mat[i, j+Ngauss+1] = V12(E, p[1][i], p[2][j]) #+ Ctct_12(g_C)
            mat[i+Ngauss+1, j+Ngauss+1] = V22(E, p[2][i], p[2][j]) #+ Ctct_22(g_C)
            # if i == Ngauss + 1 && i == j
            #     println(x0[1])
            #     println(mat[i, j])
            # end
        end
    end
    return mat
end

function tmat(Λ, E, Ngauss)
    G = gmat(Λ, E, Ngauss)
    V = vmat(Λ, E , Ngauss, -1/4)
    return inv(I - V * G) * V
end

function gaussC(N,a,b)
    xx0,ww0=gauss(N,0.0,1.0)
    return [(b-a).*xx0.+a,(b-a).*ww0]
end

function quadgaussC(f, x, w)
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    for i in eachindex(x)
        res += f(x[i]) * w[i]
    end
    return res
end

facπ = g_pi^2 / f_pi^2 / 24
zi=0.3
nn=16
nn2=24
const xxpiup=[gaussC(nn,-1,-1+zi*im)[1] ;gaussC(nn2,-1+zi*im,1+zi*im)[1] ; gaussC(nn,1+zi*im,1)[1]];
const wwpiup=[gaussC(nn,-1,-1+zi*im)[2] ;gaussC(nn2,-1+zi*im,1+zi*im)[2] ; gaussC(nn,1+zi*im,1)[2]];
const xxpidown=[gaussC(nn,-1,-1-zi*im)[1] ;gaussC(nn2,-1-zi*im,1-zi*im)[1] ; gaussC(nn,1-zi*im,1)[1]];
const wwpidown=[gaussC(nn,-1,-1-zi*im)[2] ;gaussC(nn2,-1-zi*im,1-zi*im)[2] ; gaussC(nn,1-zi*im,1)[2]];
const xxz,wwz=gauss(nn2,-1.0,1.0);
const xx0ii,ww0ii=gauss(nn2,0.0,1.0);

Epi(z,p1,p2,m0)=sqrt(p1^2+p2^2-2*p1*p2*z+m0^2)
Dij(E,z,p1,p2,mi,mj,m0)=E-(mi+p1^2/(2*(mi)))-(mj+p2^2/(2*(mj)))-Epi(z,p1,p2,m0)+im*ϵ
z0(E,m,p1,p2,m0)=((E - 2m - (p1^2 + p2^2)/(2m))^2 - m0^2 - p1^2 - p2^2)/(-2p2*p1)
z0E(p1,p2,m0)=(p1^2+p2^2+m0^2)/(2p2*p1)

function Vπu(E,p1,p2,m1,gam1,m2,gam2,m3,gam3,m4,gam4,m0,fac)  # this is only valid for real axis of 11 RS 
    D1(z)=Dij(E,z,p1,p2,m1-im*gam1/2,m3-im*gam3/2,m0)
    D2(z)=Dij(E,z,p1,p2,m2-im*gam2/2,m4-im*gam4/2,m0)
    D11int(z)=facπ*fac*(1/D1(z)+1/D2(z))/(2*Epi(z,p1,p2,m0))*(p1^2+p2^2-2*p1*p2*z)
    _z0=z0(E,m1-im*gam1/2,p1,p2,m0)
    _z0E=z0E(p1,p2,m0)
    frame=0.2
#     println([_z0,_z0E])
#     println(D11int(0.1))
    if abs(real(_z0))>1||abs(imag(_z0))>frame
        if abs(real(_z0E))>1||abs(imag(_z0E))>frame
            return quadgaussC(D11int,xxz,wwz)
        elseif imag(_z0E)>0
            return quadgaussC(D11int,xxpidown,wwpidown)
        else
            return quadgaussC(D11int,xxpiup,wwpiup)
        end
    elseif imag(_z0)>0
        if abs(real(_z0E))>1||abs(imag(_z0E))>frame
            return quadgaussC(D11int,xxpidown,wwpidown)
        elseif imag(_z0E)>0
            return quadgaussC(D11int,xxpidown,wwpidown)
        else
            _z01,_z02= real(_z0)<real(_z0E) ? [real(_z0)-0.1im,real(_z0E)+0.1im] : [real(_z0E)+0.1im,real(_z0)-0.1im]
            return (quadgaussC(D11int,(_z01+1).*xx0ii.-1,(_z01+1).*ww0ii)
                +quadgaussC(D11int,(_z02-_z01).*xx0ii.+_z01,(_z02-_z01).*ww0ii)
                +quadgaussC(D11int,(1-_z02).*xx0ii.+_z02,(1-_z02).*ww0ii))
        end
    else
        if abs(real(_z0E))>1||abs(imag(_z0E))>frame
            return quadgaussC(D11int,xxpiup,wwpiup)
        elseif imag(_z0E)>0
             _z01,_z02= real(_z0)<real(_z0E) ? [real(_z0)+0.1im,real(_z0E)-0.1im] : [real(_z0E)-0.1im,real(_z0)+0.1im]
            return (quadgaussC(D11int,(_z01+1).*xx0ii.-1,(_z01+1).*ww0ii)
                +quadgaussC(D11int,(_z02-_z01).*xx0ii.+_z01,(_z02-_z01).*ww0ii)
                +quadgaussC(D11int,(1-_z02).*xx0ii.+_z02,(1-_z02).*ww0ii))
        else
            return quadgaussC(D11int,xxpiup,wwpiup)
        end
    end
end

V11(E, p, pprime) = Vπu(E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star, gamma_B_star, m_B, 0, m_pi, 3) + Vπu(E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star, gamma_B_star, m_B, 0, m_eta, 1/3)
V12(E, p, pprime) = Vπu(E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_B_star, gamma_B_star, m_B, 0, m_K, 2^(3/2))
V21(E, p, pprime) = Vπu(E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_K, 2^(3/2))
V22(E, p, pprime) = Vπu(E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_eta, 4/3)

# V11(E, p, pprime) = Vπu(E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star, gamma_B_star, m_B, 0, m_pi, 1)
# V12(E, p, pprime) = Vπu(E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_B_star, gamma_B_star, m_B, 0, m_pi, 1)
# V21(E, p, pprime) = Vπu(E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_pi, 1)
# V22(E, p, pprime) = Vπu(E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_pi, 1)


function detImVG(Λ, E, Ngauss)
    G = gmat(Λ, E, Ngauss)
    V = vmat(Λ, E, Ngauss, -1/4)
    return det(I - V*G)
end
