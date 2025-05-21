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
    ylims!(0,1e4)
    savefig("onshellT.png")
    cfree(reinterpret(Ptr{Cvoid}, otr))
end

function cfree(ptr::Ptr{Cvoid})
    ccall(Libdl.dlsym(libscript, :Free), Cvoid, (Ptr{Cvoid},), ptr)
end

E = -2.5:0.01:0.85
# println("in julia")
# println(collect(E))
conshellT(collect(E), length(E), 200)
