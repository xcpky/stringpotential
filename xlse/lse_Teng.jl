using QuadGK, NLsolve
using CSV, DataFrames, LaTeXStrings;
# ENV["COLUMNS"] = 160;


using Plots;
# default(frame=:box, minorticks=5);
#using BenchmarkTools;

# PyPlot.matplotlib.rc("mathtext", fontset="cm")        #computer modern font
# PyPlot.matplotlib.rc("font", family="serif", size=11)  #font similar to LaTeX

# using IMinuit

# using StaticArrays
# using LinearAlgebra, StaticArrays, NonlinearSolve#,Distributions, Distances
# using Combinatorics: combinations
# using Interpolations
using LinearAlgebra

using SpecialFunctions

# include("hadronmasses.jl");
include("constants.jl")

# const ϵ=1e-9;
λ(x, y, z) = x^2 + y^2 + z^2 - 2x * y - 2y * z - 2z * x

function qsq(E, m1, m2)
    return λ(E^2, m1^2, m2^2) / (4E^2)
end

function xsqrt(x)  # the cut is along positive real aixs
    imag(x) >= 0 ? sqrt(x + 0im) : -sqrt(x - 0im)
end

const μ1 = m11 * m12 / (m11 + m12);
const Fπ = 0.092;
const facπ = g^2 / 12.0 / Fπ^2;
const Γc = 83.4e-6; #D*+ width
const Γ0 = 55.3e-6; #D*0 width
const Λ = 4.0;

# D*0 width
Gam1(E, p) = 1e-9;

Gam1(m11 + m12, 0)

# sampling by gauss
function gaussC(N, a, b)
    xx0, ww0 = gauss(N, 0.0, 1.0)
    return [(b - a) .* xx0 .+ a, (b - a) .* ww0]
end

# for integration
function quadgaussC(f, x, w)
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    for i in eachindex(x)
        res += f(x[i]) * w[i]
    end
    return res
end

# define some integrate paths
zi = 0.3
nn = 6
nn2 = 12
const xxpiup = [gaussC(nn, -1, -1 + zi * im)[1]; gaussC(nn2, -1 + zi * im, 1 + zi * im)[1]; gaussC(nn, 1 + zi * im, 1)[1]];
const wwpiup = [gaussC(nn, -1, -1 + zi * im)[2]; gaussC(nn2, -1 + zi * im, 1 + zi * im)[2]; gaussC(nn, 1 + zi * im, 1)[2]];
const xxpidown = [gaussC(nn, -1, -1 - zi * im)[1]; gaussC(nn2, -1 - zi * im, 1 - zi * im)[1]; gaussC(nn, 1 - zi * im, 1)[1]];
const wwpidown = [gaussC(nn, -1, -1 - zi * im)[2]; gaussC(nn2, -1 - zi * im, 1 - zi * im)[2]; gaussC(nn, 1 - zi * im, 1)[2]];
const xxz, wwz = gauss(nn2, -1.0, 1.0);
const xx0ii, ww0ii = gauss(nn2, 0.0, 1.0);

plot(xxpiup, label="up")
plot!(xxpidown, label="down")

# potential for one pion exchange

Epi(z, p1, p2, m0) = sqrt(p1^2 + p2^2 - 2 * p1 * p2 * z + m0^2)
Dij(E, z, p1, p2, mi, mj, m0) = E - (mi + p1^2 / (2 * (mi))) - (mj + p2^2 / (2 * (mj))) - Epi(z, p1, p2, m0) + im * ϵ
z0(E, m, p1, p2, m0) = ((E - 2m - (p1^2 + p2^2) / (2m))^2 - m0^2 - p1^2 - p2^2) / (-2p2 * p1)
z0E(p1, p2, m0) = (p1^2 + p2^2 + m0^2) / (2p2 * p1)

function Vπu(E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0, fac)  # this is only valid for real axis of 11 RS 
    D1(z) = Dij(E, z, p1, p2, m1 - im * gam1 / 2, m3 - im * gam3 / 2, m0)
    D2(z) = Dij(E, z, p1, p2, m2 - im * gam2 / 2, m4 - im * gam4 / 2, m0)
    D11int(z) = facπ * fac * (1 / D1(z) + 1 / D2(z)) / (2 * Epi(z, p1, p2, m0)) * (p1^2 + p2^2 - 2 * p1 * p2 * z)
    _z0 = z0(E, m1 - im * gam1 / 2, p1, p2, m0)
    _z0E = z0E(p1, p2, m0)
    frame = 0.2
    #     println([_z0,_z0E])
    #     println(D11int(0.1))
    if abs(real(_z0)) > 1 || abs(imag(_z0)) > frame
        if abs(real(_z0E)) > 1 || abs(imag(_z0E)) > frame
            return quadgaussC(D11int, xxz, wwz)
        elseif imag(_z0E) > 0
            return quadgaussC(D11int, xxpidown, wwpidown)
        else
            return quadgaussC(D11int, xxpiup, wwpiup)
        end
    elseif imag(_z0) > 0
        if abs(real(_z0E)) > 1 || abs(imag(_z0E)) > frame
            return quadgaussC(D11int, xxpidown, wwpidown)
        elseif imag(_z0E) > 0
            return quadgaussC(D11int, xxpidown, wwpidown)
        else
            _z01, _z02 = real(_z0) < real(_z0E) ? [real(_z0) - 0.1im, real(_z0E) + 0.1im] : [real(_z0E) + 0.1im, real(_z0) - 0.1im]
            return (quadgaussC(D11int, (_z01 + 1) .* xx0ii .- 1, (_z01 + 1) .* ww0ii)
                    + quadgaussC(D11int, (_z02 - _z01) .* xx0ii .+ _z01, (_z02 - _z01) .* ww0ii)
                    + quadgaussC(D11int, (1 - _z02) .* xx0ii .+ _z02, (1 - _z02) .* ww0ii))
        end
    else
        if abs(real(_z0E)) > 1 || abs(imag(_z0E)) > frame
            return quadgaussC(D11int, xxpiup, wwpiup)
        elseif imag(_z0E) > 0
            _z01, _z02 = real(_z0) < real(_z0E) ? [real(_z0) + 0.1im, real(_z0E) - 0.1im] : [real(_z0E) - 0.1im, real(_z0) + 0.1im]
            return (quadgaussC(D11int, (_z01 + 1) .* xx0ii .- 1, (_z01 + 1) .* ww0ii)
                    + quadgaussC(D11int, (_z02 - _z01) .* xx0ii .+ _z01, (_z02 - _z01) .* ww0ii)
                    + quadgaussC(D11int, (1 - _z02) .* xx0ii .+ _z02, (1 - _z02) .* ww0ii))
        else
            return quadgaussC(D11int, xxpiup, wwpiup)
        end
    end
end

# sampling of discreet momentun in [0,Lambda]

xx110, ww110 = gaussC(30, 0, 0.2)
xx120, ww120 = gaussC(10, 0.2, Λ)
xx1 = [xx110; xx120]
ww1 = [ww110; ww120]
xx1, ww1 = gauss(64, 0, Λ)
xx1 = convert(Array{Complex}, xx1)
n1 = length(xx1)

vmat = complex(zeros(n1 + 1, n1 + 1));
gmat = complex(zeros(n1 + 1, n1 + 1));
tmat = complex(zeros(n1 + 1, n1 + 1));

using JLD

function Vmat!(E, rs, c0x, c1x, vmat)
    d0sw = 0

    # compute on-shell momentum
    dE1 = E - m11 - m12
    q01 = xsqrt(2μ1 * dE1)
    x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, q01) / 2)) #irratate once
    x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, x01) / 2)) #irratate twice
    x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, x01) / 2)) #irratate 3 times
    x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, x01) / 2)) #irratate 4 times

    if rs == 1
        nothing
    else
        x01 = -x01
    end

    xx11 = [xx1; x01]
    nn1 = n1 + 1

    for i = 1:nn1
        for j = 1:nn1
            # vmat[i, j] = (c0x + c1x) / 2 + Vπu(E, xx11[i], xx11[j], m11, d0sw, m12, Gam1(E, xx11[i]), m11, d0sw, m12, Gam1(E, xx11[j]), m50, 0.5)
			vmat[i, j] = xx11[i]*xx11[j]
        end
    end


    return vmat

end

function Gmat!(E, rs, gam1, gmat)

    # compute on-shell momentum
    dE1 = E - m11 - m12
    # q01 = xsqrt(2μ1 * dE1)
    # x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, q01) / 2)) #irratate once
    # x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, x01) / 2)) #irratate twice
    # x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, x01) / 2)) #irratate 3 times
    # x01 = xsqrt(2μ1 * (dE1 + 1im * Gam1(E, x01) / 2)) #irratate 4 times
    x01 = xsqrt(2mu[1] * dE1)

    if rs == 1
        nothing
    else
        x01 = -x01
    end

    xx11 = [xx1; x01]

    integrand1(q) = 1 / (dE1 - q^2 / 2μ1 + 1im * Gam1(E, x01) / 2 + 1im * gam1 / 2)

    for i = 1:n1
        gmat[i, i] = 1 / (2 * π^2) * xx1[i]^2 * ww1[i] / (dE1 - xx1[i]^2 / 2μ1 + 1im * Gam1(E, xx1[i]) / 2 + 1im * gam1 / 2)
    end
    gmat[n1+1, n1+1] = 1 / (2 * π^2) * x01 * (μ1 * (-im * π + log((Λ + x01) / (Λ - x01))) - x01 * quadgaussC(integrand1, xx1, ww1))

    return gmat
end

function Tmat!(E, rs, c0x, c1x, gam1, vmat, gmat, tmat)
    Vmat!(E, rs, c0x, c1x, vmat)
    Gmat!(E, rs, gam1, gmat)
    tmat = inv(I - vmat * gmat) * vmat

    # return det(I - vmat * gmat)
    return tmat[end, end]
end

Erange = LinRange(m_Xb11P - 0.1, m_Xb14P + 0.1, 1000) .+ (m11 + m12)
function ot()
    t = [Tmat!(Erange[i], 1, 0, 0, 0, vmat, gmat, tmat) for i in eachindex(Erange)]
    plot(dpi=400)
    yup = 1.2maximum(abs.(t))
    ylw = 0.8minimum(abs.(t))

    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    vline!([m_Xb14P], s=:dash, label=L"\chi_{b1}(4P)")
    plot!(Erange .- (m11 + m12), abs.(t))
	xlabel!("E/GeV")
    ylims!(ylw, yup)
	xlims!(Erange[1]-(m11+m12), Erange[end]-(m11+m12))
    savefig("onshellTeng.png")
end

function trg()
    g = [tr(Gmat!(Erange[i], 1, 0, gmat)) for i in eachindex(Erange)]
    plot(dpi=400)
    vline!([0, m21 + m22 - m11 - m12], label="thresholds", s=:dash)
    plot!(Erange .- (m11 + m12), real.(g), label="real")
    plot!(Erange .- (m11 + m12), imag.(g), label="imag")
    ylims!(-1.5, 0.5)
    savefig("trgeng.png")
end

# Gmat!(0.1 + m11 + m12, 1, 0, gmat);
trg()
ot()
