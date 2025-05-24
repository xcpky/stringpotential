using Libdl
# Load the shared library
const libscript = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/libscript.so"))
using Plots
using LaTeXStrings

function conshellT(E::Vector{Cdouble}, len, pNgauss)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint), E, len, pNgauss)
    ot = copy(transpose(unsafe_wrap(Array, otr, (len, 4), own=false)))
    plot(E, abs.(ot[1,:]), label=L"$T_{11}$", dpi=300)
    plot!(E, abs.(ot[2,:]), label=L"$T_{12}$")
    plot!(E, abs.(ot[4,:]), label=L"$T_{22}$")
    # ylims!(0,1e4)
    savefig("onshellT.png")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function detImVG(E::Vector{Cdouble}, len, pNgauss)
    @time dtr = ccall(Libdl.dlsym(libscript, :Det), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint), E, len, pNgauss)
    Det = copy(unsafe_wrap(Array, dtr, len, own=false))
    plot(E, real.(Det), label=L"det($1-VG$)",dpi=300)
    plot!(E, imag(Det), label=L"im")
    savefig("det.png")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return Det
end

function cfree(ptr::Ptr{Cvoid})
    ccall(Libdl.dlsym(libscript, :Free), Cvoid, (Ptr{Cvoid},), ptr)
end

# E = 0.499:0.000013:0.5035
E = 0.89:0.000001:0.91
# println("in julia")
# println(collect(E))
ot = conshellT(collect(E), length(E), 64)
# Det = detImVG(collect(E), length(E), 200)
