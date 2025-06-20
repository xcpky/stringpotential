using Libdl
# Load the shared library
const libscript = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/libscript.so"))
using Plots
using LaTeXStrings

function conshellT(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    ot = copy(transpose(unsafe_wrap(Array, otr, (len, 4), own=false)))
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    # ylims!(0,1e4)
    savefig("onshellT.png")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function detImVG(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Det), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    Det = copy(unsafe_wrap(Array, dtr, len, own=false))
    # plot(E, real.(Det), label=L"det($1-VG$)",dpi=300)
    plot(E, abs.(Det), label=L"|det($1-VG$)|", dpi=300)
    # ylims!(0,4)
    savefig("det.png")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return Det
end

function Both(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Both), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    data = copy((unsafe_wrap(Array, dtr, (5, len), own = false)))
    ot = data[2:5, :]
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    # ylims!(0,1e4)
    # savefig("onshellT.png")
    plot!(E, abs.(data[1, :] .- 1), label=L"|det($1-VG$)|", dpi=300)
    # ylims!(0,1.4)
    ylims!(0, 1.2)
    level = getEvec(10, Lambda, epsilon)
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline!(vls, s=:dash, c=:grey, label = "quark model")
    # ylims!(0, 2.1)
    savefig("det.png")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
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

# E = 1.48:0.00001:1.499
# E = LinRange(1.37, 1.3725, 1000)
E = LinRange(-2.7325, 0.732, 5000)
# E = 0.399999999:0.00000000001:0.400000001
# E = 0:0.01:1.8
# println("in julia")
# println(collect(E))
ot = conshellT(collect(E), length(E), 64, 4, 1e-6)
# Det = detImVG(collect(E), length(E), 64, 4, 1e-6)
# data = Both(collect(E), length(E), 64, 4, 1e-6)
# E = getEvec(64, 4, 1e-6)
