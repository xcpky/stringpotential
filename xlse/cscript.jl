using Libdl
# Load the shared library
const libscript = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/libscript.so"))
using Plots
using LaTeXStrings

include("constants.jl")

function conshellT(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    ot = copy(transpose(unsafe_wrap(Array, otr, (len, 4), own=false)))
    upper = min(1e4, maximum(abs.(ot)))
    level = getEvec(10, Lambda, epsilon)
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # ylims!(0, upper)
    # ylims!(0,10)
    savefig("onshellT.png")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function conshellG(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellG), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    ot = copy(transpose(unsafe_wrap(Array, otr, (len, 4), own=false)))
    upper = min(1e4, maximum(abs.(ot)))
    level = getEvec(10, Lambda, epsilon)
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    ylims!(0, upper)
    # ylims!(0,10)
    savefig("onshellG.png")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function conshellV(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellV), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    ot = copy(transpose(unsafe_wrap(Array, otr, (len, 4), own=false)))
    upper = min(1e4, maximum(abs.(ot)))
    level = getEvec(10, Lambda, epsilon)
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # ylims!(0, upper)
    # ylims!(0,10)
    savefig("onshellV.png")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function detImVG(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Det), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    Det = copy(unsafe_wrap(Array, dtr, len, own=false))
    # plot(E, real.(Det), label=L"det($1-VG$)",dpi=300)
    plot(E, abs.(Det), label=L"|det($1-VG$)|", dpi=300)
    savefig("det.png")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return Det
end

function Both(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Both), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    data = copy((unsafe_wrap(Array, dtr, (5, len), own=false)))
    ot = data[2:5, :]
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    # ylims!(0,1e4)
    # savefig("onshellT.png")
    plot!(E, abs.(data[1, :]), label=L"|det($1-VG$)|", dpi=300)
    # ylims!(0,1.4)
    # ylims!(0, 1.2)
    level = getEvec(10, Lambda, epsilon)
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    # ylims!(0, 2.1)
    savefig("det.png")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return data
end

function Poles(Er::Vector{Cdouble}, rlen, Ei::Vector{Cdouble}, ilen, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Poles), Ptr{ComplexF64}, (Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Csize_t, Csize_t, Cdouble, Cdouble), Er, rlen, Ei, ilen, pNgauss, Lambda, epsilon)
    data = copy(unsafe_wrap(Array, dtr, rlen * ilen * 4, own=false))
    data = transpose(reshape(data, (rlen * ilen, 4)))
    poles = [filter(!isnan, data[i, :]) for i in 1:4]
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    scatter(real.(poles[1]), imag.(poles[1]), dpi=300, label=L"RS_{--}", markersize=3, markerstrokewidth=0,legend=:outertopright)
    scatter!(real.(poles[2]), imag.(poles[2]), dpi=300, label=L"RS_{+-}", markersize=3, markerstrokewidth=0)
    scatter!(real.(poles[3]), imag.(poles[3]), dpi=300, label=L"RS_{-+}", markersize=3, markerstrokewidth=0)
    scatter!(real.(poles[4]), imag.(poles[4]), dpi=300, label=L"RS_{++}", markersize=3, markerstrokewidth=0)
    level = getEvec(10, Lambda, epsilon)
    upper = -1e6
    lower = 1e6
    for i in 1:4
        upper = max(maximum(real.(poles[i])), upper)
        lower = min(minimum(real.(poles[i])), lower)
    end
    # println(level)
    vls = filter(e -> e > lower - 0.2 && e < upper + 0.2, level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    xlabel!("Re")
    ylabel!("Im")
    savefig("pole.png")
    return data
end

function getEvec(pNgauss, Lambda, epsilon)
    etr = ccall(Libdl.dlsym(libscript, :Evec), Ptr{Cdouble}, (Cuint, Cdouble, Cdouble), pNgauss, Lambda, epsilon)
    E = copy(unsafe_wrap(Array, etr, 64, own=false))
    return E
end

function cfree(ptr::Ptr{Cvoid})
    ccall(Libdl.dlsym(libscript, :Free), Cvoid, (Ptr{Cvoid},), ptr)
end

epsi = 1e-5
Lambda = 4
pNgauss = 128
data = Nothing
if "--poles" in ARGS
    Er = LinRange(-0.8, 0.5, 64)
    Ei = LinRange(-0.01, 0.01, 16)
    data = Poles(collect(Er), length(Er), collect(Ei), length(Ei), 64, Lambda, epsi)
end
onshellRange = LinRange(-1, -0.76, 5000)

if "--onshellT" in ARGS
    # E = 1.48:0.00001:1.499
    E = onshellRange
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    # println("in julia")
    # println(collect(E))
    ot = conshellT(collect(E), length(E), pNgauss, Lambda, epsi)
    # og = conshellG(collect(E), length(E), 100, 4, epsi)
    # Det = detImVG(collect(E), length(E), 64, 4, epsi)
    # data = Both(collect(E), length(E), 64, 4, epsi)
    # E = getEvec(64, 4, 1e-6)
end
if "--onshellG" in ARGS
    # E = 1.48:0.00001:1.499
    E = onshellRange
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    # println("in julia")
    # println(collect(E))
    # ot = conshellG(collect(E), length(E), 100, 4, epsi)
    og = conshellG(collect(E), length(E), pNgauss, Lambda, epsi)
    # Det = detImVG(collect(E), length(E), 64, 4, epsi)
    # data = Both(collect(E), length(E), 64, 4, epsi)
    # E = getEvec(64, 4, 1e-6)
end
if "--onshellV" in ARGS
    # E = 1.48:0.00001:1.499
    E = onshellRange
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    ov = conshellV(collect(E), length(E), pNgauss, Lambda, epsi)
end
