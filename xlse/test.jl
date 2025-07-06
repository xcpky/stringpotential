using Libdl

# Load the shared library
const libscript = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/libscript.so"))
pNgauss = 40
Lambda = 4
epsilon = 1e-9
E = 0.5
n = 2*(pNgauss + 1)
lse = ccall(Libdl.dlsym(libscript, :lse_malloc), Ptr{Cvoid}, (Csize_t, Cdouble, Cdouble), pNgauss, Lambda, epsilon)
ccall(Libdl.dlsym(libscript, :lse_refresh), Cvoid, (Ptr{Cvoid}, ComplexF64, Ptr{Cdouble}, Cuint), lse, E, [0., 0, 0, 0], 3)
ccall(Libdl.dlsym(libscript, :lse_vmat), Cint, (Ptr{Cvoid},), lse)
ptr = ccall(Libdl.dlsym(libscript, :lse_get_v_data), Ptr{ComplexF64}, (Ptr{Cvoid},), lse)
vv = transpose(copy(unsafe_wrap(Array, ptr, (n, n), own=false)))
