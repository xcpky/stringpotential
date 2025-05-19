using Plots: Legend
# Load the shared library
include("lse.jl")
using Plots
using Base.Threads
using .CLSE
using QuadGK
using LinearAlgebra
using StatsPlots
using Measures
using StatsBase
using Serialization
using LaTeXStrings
Ngauss = 200
eps = 1e-7
E = -2.5:0.01:0.85
Λ = 4
xi, wi = gauss(Ngauss, 0, Λ)
function trOnshellT(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_compute(lse, E, rs)
    T = lse_get_t_data(lse)
    return tr(T[1:Ngauss+1, 1:Ngauss+1]), tr(T[1:Ngauss+1, Ngauss+2:end]), tr(T[Ngauss+2:end, 1:Ngauss+1]), tr(T[Ngauss+2:end, Ngauss+2:end])
end

function OnshellT(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_compute(lse, E, rs)
    T = lse_get_t_data(lse)
    return T[Ngauss+1, Ngauss+1], T[Ngauss+1, end], T[end, Ngauss+1], T[end, end]
end

function OnshellG(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_refresh(lse, E, rs)
    lse_gmat(lse)
    G = lse_get_g_data(lse)
    return [tr(G[1:Ngauss+1, 1:Ngauss+1]), tr(G[Ngauss+2:2*Ngauss+2, Ngauss+2:2*Ngauss+2])]
end
function onshell(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_refresh(lse, E, 1)
    lse_gmat(lse)
    G = lse_get_g_data(lse)
    return G[Ngauss+1, Ngauss+1], G[Ngauss+1, end], G[end, Ngauss+1], G[end]
end
function OnshellV(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_refresh(lse, E, rs)
    lse_vmat(lse)
    V = lse_get_v_data(lse)
    return V[Ngauss+1, Ngauss+1], V[Ngauss+1, end], V[end, Ngauss+1], V[end]
end
function OnshellPsi(lse::Ptr{LSE}, E::ComplexF64, n::Int, rs::Int)
    lse_refresh(lse, E, rs)
    psi = lse_get_psi(lse)
    return psi[n, end-2], psi[n, end-1]
end

function onshelliIVG(lse::Ptr{LSE}, E::ComplexF64)
    lse_compute(lse, E, 1)
    ivg = lse_get_ivg_data(lse)
    iivg = inv(ivg)
    return iivg[Ngauss+1, Ngauss+1], iivg[Ngauss+1, end], iivg[end, Ngauss+1], iivg[end]
end

function testiivg(lse::Ptr{LSE}, E::ComplexF64)
    lse_compute(lse, E, 1)
    ivg = lse_get_ivg_data(lse)
    return inv(ivg)
end

function detiIVG(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_compute(lse, E, rs)
    ivg = lse_get_ivg_data(lse)
    return det(inv(ivg))
end

function detG(lse::Ptr{LSE}, E)
    lse_refresh(lse, E, 1)
    lse_gmat(lse)
    return det(copy(lse_get_g_data(lse)[1:Ngauss+1, 1:Ngauss+1]))
end

function detV(lse::Ptr{LSE}, E)
    lse_refresh(lse, E, 1)
    lse_vmat(lse)
    return det(copy(lse_get_v_data(lse)))
end

function detIVG(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_compute(lse, E, rs)
    return det(lse_get_ivg_data(lse))
end

function detT(lse::Ptr{LSE}, E::ComplexF64, rs::Int)
    lse_compute(lse, E, rs)
    t = lse_get_t_data(lse)
    return abs(det(t))
end

function detImVG!(lse::Ptr{LSE}, F, x)
    w = x[1] + im * x[2]
    det = lse_detImVG(lse, w)
    F[1] = real(det)
    F[2] = imag(det)
end

using NLsolve
function poles(lse::Ptr{LSE}, init_val=[0, 0])
    sol = nlsolve((F, x) -> detImVG!(lse, F, x), init_val)
    # println(sol.residual_norm)
    return sol.residual_norm < 1e-5 ? Complex(sol.zero[1], sol.zero[2]) : NaN + im * NaN
end

function lse_iivg(lse::Ptr{LSE}, E::Real)
    lse_refresh(lse, E + 0im, 1)
    return copy(lse_get_iivg_data(lse))
end

# Example usage
# wf = wf_new(0, 20, 40)
# c_matrix = get_c_solution_matrix(wf)
# e_vector = get_e_solution_vector(wf)
# println("C solution matrix: ", c_matrix)
lse = lse_malloc(Ngauss, Λ, eps)
@time lse_compute(lse, -0.2 + 0im, 1)
# E = -0.578802975:0.0000000001:-0.57880295
# E = -3:0.0003:-2

if "--testV" in ARGS
    v = copy(lse_get_v_data(lse))
end

if "--testWF" in ARGS
    wf = wf_new(1, 20, 64)
    r = collect(range(0, 7, 1000))
    psi = psi_n_batch(wf, xi, UInt(4), 0.0)
    plot(r, imag.(psi_n_ft_batch(wf, r, UInt(1))), dpi=300)
    plot!(r, imag.(psi_n_ft_batch(wf, r, UInt(2))), dpi=300)
    plot!(r, imag.(psi_n_ft_batch(wf, r, UInt(3))), dpi=300)
    plot!(r, imag.(psi_n_ft_batch(wf, r, UInt(4))), dpi=300)
    plot!(r, imag.(psi_n_ft_batch(wf, r, UInt(5))), dpi=300)
    savefig("psi.png")
end

if "--testPsi" in ARGS
    lse_refresh(lse, 0.2 + 0im, 1)
    data = read("psi.dat")
    psi = transpose(reshape(reinterpret(ComplexF64, data), (Ngauss+2, 30)))
    plot(xi, real.(psi[1, 1:Ngauss]), dpi=300)
    plot!(xi, real.(psi[2, 1:Ngauss]), dpi=300)
    plot!(xi, real.(psi[3, 1:Ngauss]), dpi=300)
    plot!(xi, real.(psi[4, 1:Ngauss]), dpi=300)
    # xlims!(1, 5)
    savefig("psi.png")
    # for i in 1:30
    #     local quad = 0
    #     for j in 1:Ngauss
    #         quad += abs(psi[i, j]) * wi[j]
    #     end
    #     println(quad)
    # end
    # println(psi)
end

if "--onshellPsi" in ARGS
    op = OnshellPsi.(lse, Complex.(E), 1, 1)
    psi0 = [op[i][1] for i in eachindex(E)]
    psi1 = [op[i][2] for i in eachindex(E)]
end

if "--onshellG" in ARGS
    tce = OnshellG.(lse, Complex.(E), 1)
    g11 = [tce[i][1] for i in eachindex(E)]
    g22 = [tce[i][2] for i in eachindex(E)]
    # plot(E, real.(g11), dpi=300)
    # plot!(E, imag.(g11))
    # plot!(E, real.(g22))
    # plot!(E, imag.(g22))
    plot(E, abs.(g11), dpi=300)
    plot!(E, abs.(g22))
    savefig("onshellG.png")
end

if "--onshellV" in ARGS
    oV = OnshellV.(lse, Complex.(E), 1)
    oV11 = [abs(oV[i][1]) for i in eachindex(E)]
    # oV12 = [abs(oV[i][2]) for i in eachindex(E)]
    # oV21 = [abs(oV[i][3]) for i in eachindex(E)]
    oV22 = [abs(oV[i][4]) for i in eachindex(E)]
    # plot(E, oV11, dpi=300, label="abs(V_OME), channel=1")
    plot(E, oV22, dpi=300, label="abs(V_OME), channel=2")
    # ylims!(0, 1e5)
    # plot!(E, oV12, label="alpha 1 beta 2")
    # plot!(E, oV21, label="alpha 2 beta 1")
    # plot!(E, oV22, label="alpha 2 beta 2")
    savefig("onshellV.png")
end

if "--onshellT" in ARGS
    oT = Array{NTuple{4,ComplexF64}}(undef, length(E))
    if "--cached" in ARGS
        oT = deserialize("onshellT.dat")
    else
        thread_lse = Vector{Any}(undef, nthreads())
        @threads for i in 1:nthreads()
            thread_lse[i] = lse_malloc(Ngauss, 4, eps)
        end    #
        println(nthreads())
        oT = Array{NTuple{4,ComplexF64}}(undef, length(E))
        @threads for i in eachindex(E)
            lse = thread_lse[threadid()]
            oT[i] = OnshellT(lse, Complex(E[i]), 1)
        end
        serialize("onshellT.dat", oT)
    end
    # oT = OnshellT.(lse, Complex.(E), 1)
    oT11 = [abs(oT[i][1]) for i in eachindex(E)]
    oT12 = [abs(oT[i][2]) for i in eachindex(E)]
    oT21 = [abs(oT[i][3]) for i in eachindex(E)]
    oT22 = [abs(oT[i][4]) for i in eachindex(E)]
    h = max(maximum(oT11), maximum(oT12), maximum(oT21), maximum(oT22))
    plot(E, oT11, dpi=300, label=L"$\alpha=1,\beta=1$")
    # xlims!(-3, -2)
    # ylims!(0, 1e5)
    # ylims!(0, 0.75*h)
    plot!(E, oT12, label=L"$\alpha=1,\beta=2$")
    # plot!(E, oT21, label=L"$\alpha=2,\beta=1$")
    plot!(E, oT22, label=L"$\alpha=2,\beta=2$")
    xlabel!("E/GeV")
    ylabel!(L"$|T_{\alpha\beta}|GeV^{-2}$")
    savefig("onshellT.png")
end

if "--detiIVG" in ARGS
    detiivg1 = detiIVG.(lse, Complex.(E), 1)
    detiivg2 = detiIVG.(lse, Complex.(E), -1)
end

if "--detIVG" in ARGS
    # detivg1 = detIVG.(lse, Complex.(E), 1 )
    detivg1 = lse_detImVG.(lse, Complex.(E))
    # detivg2 = detIVG.(lse, E, -1, 0, 0)
    plot(E, inv.(abs.(detivg1)), dpi=300, label="inverse(det(I-VG))")
    ylims!(0, 10)
    # plot!(E, imag.(detivg1), label="imag")
    # plot!(E, real.(detivg2), label="rs- real")
    # plot!(E, imag.(detivg2), label="rs- imag")
    savefig("detIVG.png")
end

if "--detT" in ARGS
    detv = detT.(lse, Complex.(E), 1)
    plot(E, detv, dpi=300)
    savefig("detT.png")
end

if "--poles" in ARGS
    re = range(-1.5, 0, 16)
    ie = range(-1.5, 0, 16)

    # Initialize the result array first
    pole = Array{Any}(undef, length(re) * length(ie))
    thread_lse = Vector{Any}(undef, nthreads())

    @threads for i in 1:nthreads()
        thread_lse[i] = lse_malloc(Ngauss, 4, eps)
    end    #
    @threads for i in 1:length(re)*length(ie)
        # Get thread-local lse
        lse = thread_lse[threadid()]

        # Convert linear index to 2D indices
        e_idx = div(i - 1, length(ie)) + 1  # Fixed indexing calculation
        f_idx = mod1(i, length(ie))

        # Get the actual elements
        e = re[e_idx]
        f = ie[f_idx]

        # Compute the result and store it
        tmp = poles(lse, [e, f])
        # println(typeof(pole[i]))
        pole[i] = abs(tmp) < 10 && real(tmp) < 0.2 && abs(imag(tmp)) < 0.5 ? tmp : NaN + im * NaN
    end

    # Free all thread-local lse objects
    for lse in thread_lse
        lse_free(lse)
    end
    # pole = [poles([e, 0]) for e in E]
    # filtered_pole = filter(z -> abs(z) <= 10, pole)
    # filtered_pole = filter(z -> real(z) <= 1.1, filtered_pole)

    # using Plots
    # scatter(real.(pole), imag.(pole), dpi=300,
    #     label="det(I-VG)=0",
    #     markersize=3,
    #     markerstrokewidth=0)
    # xlabel!("Re")
    # ylabel!("Im")
    # # ylims!(-1, 1)
    # savefig("poles.png")
end

@inline function xsqrt(x)
    imag(x) >= 0 ? sqrt(x + 0im) : -sqrt(x - 0im)
end

if "--Os" in ARGS
    include("constants.jl")
    Os = Array{ComplexF64}(undef, 4, length(E))
    for e in eachindex(E)
        dE = -(delta .- E[e])
        p0 = xsqrt.(2 .* mu .* dE)
        Os[1, e] = O_00(E[e] + 0im, p0[1], p0[2], m_eta)
        Os[2, e] = O_01(E[e] + 0im, p0[1], p0[2], m_eta)
        Os[3, e] = O_10(E[e] + 0im, p0[1], p0[2], m_eta)
        Os[4, e] = O_11(E[e] + 0im, p0[1], p0[2], m_eta)
    end
    plot(E, abs.(Os[1, :]), label="absolute value O_00", dpi=300)
    savefig("O_00.png")
    plot(E, abs.(Os[2, :]), label="absolute value O_01", dpi=300)
    savefig("O_01.png")
    plot(E, abs.(Os[3, :]), label="absolute value O_10", dpi=300)
    savefig("O_10.png")
    plot(E, abs.(Os[4, :]), label="absolute value O_11", dpi=300)
    savefig("O_11.png")
end

if "--onshelliIVG" in ARGS
    iivg = onshelliIVG.(lse, Complex.(E))
    iivg11 = [abs(iivg[i][1]) for i in eachindex(E)]
    iivg12 = [abs(iivg[i][2]) for i in eachindex(E)]
    iivg21 = [abs(iivg[i][3]) for i in eachindex(E)]
    iivg22 = [abs(iivg[i][4]) for i in eachindex(E)]
    plot(E, abs.(iivg11), label="onshell inverse(I-VG) alpha 1 beta 1", dpi=300)
    plot!(E, abs.(iivg12), label="onshell inverse(I-VG) alpha 1 beta 2")
    plot!(E, abs.(iivg21), label="onshell inverse(I-VG) alpha 2 beta 1")
    plot!(E, abs.(iivg22), label="onshell inverse(I-VG) alpha 2 beta 2")
    savefig("onshelliIVG.png")
end

if "--onshellGV" in ARGS
    oV = OnshellV.(lse, Complex.(E), 1)
    oG = onshell.(lse, Complex.(E), 1)
    ov11 = [abs(oV[i][1]) for i in eachindex(E)]
    ov22 = [abs(oV[i][4]) for i in eachindex(E)]
    og11 = [abs(oG[i][1]) for i in eachindex(E)]
    og22 = [abs(oG[i][4]) for i in eachindex(E)]
    # plot(E, ov11, dpi=300, label="abs(V_OME), channel=2")
    plot(E, og11, dpi=300, label="abs(G), channel=2")
    savefig("onshellGV.png")
end


if "--detGV" in ARGS
    detv = detV.(lse, Complex.(E))
    detg = detG.(lse, Complex.(E))
    plot(E, abs.(detv), dpi=300, label="det(V)")
    plot!(E, abs.(detg), dpi=300, label="det(G)")
    savefig("detGV.png")
end

if "--iIVGdist" in ARGS
    E = [-3, -0.2, -2, -1.5]
    # E = -0.2
    local i = 1
    # p = plot(layout=(3, 1), size=(800,800), dpi=300)
    for e in E
        mat = abs.(testiivg(lse, e + 0im))
        heatmap(mat,
            aspect_ratio=:equal,  # Keep cells square
            label="E=$(e)",
            dpi=300,
            c=:coolwarm)  # Color scheme    savefig("invI_VG_dist.png")
        title!("abs(I-VG) E=$(e)")
        i += 1

        savefig("heatmap$(e).png")
    end
    # savefig("heatmap.png")
end

