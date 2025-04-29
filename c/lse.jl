module CLSE

using Libdl

# Load the shared library
const liblse = Libdl.dlopen("/home/choros/code/stringpotential/c/liblse.so")

# Define struct to match the C structure
mutable struct LSE
    Ngauss::Csize_t
    Lambda::Cdouble
    epsilon::Cdouble
    E::Cdouble
    T::Ptr{Cvoid}
    V::Ptr{Cvoid}
    G::Ptr{Cvoid}
    xi::Ptr{Cdouble}
    wi::Ptr{Cdouble}
    x0::NTuple{2, Cdouble}
    table::Ptr{Cvoid}
    wf::Ptr{Cvoid}
    psi_n_mat::Ptr{Ptr{ComplexF64}}
    E_vec::Ptr{Cdouble}
end

# Function bindings
function lse_malloc(Ngauss::Integer, Lambda::Real, epsilon::Real)
    ccall(
        Libdl.dlsym(liblse, :lse_malloc),
        Ptr{LSE},
        (Csize_t, Cdouble, Cdouble),
        Ngauss, Lambda, epsilon
    )
end

function lse_compute(lse_ptr::Ptr{LSE}, E::Real)
    ccall(
        Libdl.dlsym(liblse, :lse_compute),
        Cint,
        (Ptr{LSE}, Cdouble),
        lse_ptr, E
    )
end

function lse_free(lse_ptr::Ptr{LSE})
    ccall(
        Libdl.dlsym(liblse, :lse_free),
        Cvoid,
        (Ptr{LSE},),
        lse_ptr
    )
end

function lse_refresh(lse_ptr::Ptr{LSE}, E::Real)
    ccall(
        Libdl.dlsym(liblse, :lse_refresh),
        Cvoid,
        (Ptr{LSE}, Cdouble),
        lse_ptr, E
    )
end

# Functions to get matrix data
function lse_get_t_data(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_t_data),
        Ptr{ComplexF64},
        (Ptr{LSE},),
        lse_ptr
    )
    
    rows = Ref{Cuint}(0)
    cols = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_t_size),
        Cvoid,
        (Ptr{LSE}, Ptr{Cuint}, Ptr{Cuint}),
        lse_ptr, rows, cols
    )
    
    # Convert to Julia array
    return unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false)
end

function lse_get_g_data(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_g_data),
        Ptr{ComplexF64},
        (Ptr{LSE},),
        lse_ptr
    )
    
    rows = Ref{Cuint}(0)
    cols = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_g_size),
        Cvoid,
        (Ptr{LSE}, Ptr{Cuint}, Ptr{Cuint}),
        lse_ptr, rows, cols
    )
    
    # Convert to Julia array
    return unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false)
end

function lse_gmat(lse_ptr::Ptr{LSE})
    ccall(
        Libdl.dlsym(liblse, :lse_gmat),
        Cint,
        (Ptr{LSE},),
        lse_ptr
    )
end

function lse_vmat(lse_ptr::Ptr{LSE})
    ccall(
        Libdl.dlsym(liblse, :lse_vmat),
        Cint,
        (Ptr{LSE},),
        lse_ptr
    )
end

function lse_tmat(lse_ptr::Ptr{LSE})
    ccall(
        Libdl.dlsym(liblse, :lse_tmat),
        Cint,
        (Ptr{LSE},),
        lse_ptr
    )
end

function lse_get_v_data(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_v_data),
        Ptr{ComplexF64},
        (Ptr{LSE},),
        lse_ptr
    )
    
    rows = Ref{Cuint}(0)
    cols = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_v_size),
        Cvoid,
        (Ptr{LSE}, Ptr{Cuint}, Ptr{Cuint}),
        lse_ptr, rows, cols
    )
    
    # Convert to Julia array
    return unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false)
end

export LSE, lse_free, lse_malloc, lse_compute, lse_get_t_data, lse_get_g_data, lse_get_v_data, lse_refresh, lse_gmat, lse_vmat, lse_tmat

end # module
