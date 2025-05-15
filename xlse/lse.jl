module CLSE

using Libdl

# Load the shared library
const liblse = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/liblse.so"))

# Define struct to match the C structure
mutable struct LSE
    Ngauss::Csize_t
    Lambda::Cdouble
    epsilon::Cdouble
    E::ComplexF64
    T::Ptr{Cvoid}
    V::Ptr{Cvoid}
    G::Ptr{Cvoid}
    iIVG::Ptr{Cvoid}
    xi::Ptr{Cdouble}
    wi::Ptr{Cdouble}
    x0::NTuple{2, Cdouble}
    table::Ptr{Cvoid}
    wf::Ptr{Cvoid}
    psi_n_mat::Ptr{Ptr{ComplexF64}}
    psi_raw::Ptr{ComplexF64}
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

function lse_compute(lse_ptr::Ptr{LSE}, E::ComplexF64, rs::Int)
    ccall(
        Libdl.dlsym(liblse, :lse_compute),
        Cint,
        (Ptr{LSE}, ComplexF64, Int64),
        lse_ptr, E, rs
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

function lse_refresh(lse_ptr::Ptr{LSE}, E::ComplexF64, rs::Int)
    ccall(
        Libdl.dlsym(liblse, :lse_refresh),
        Cvoid,
        (Ptr{LSE}, ComplexF64, Int64),
        lse_ptr, E, rs
    )
end

function lse_get_iivg_data(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_iivg_data),
        Ptr{ComplexF64},
        (Ptr{LSE},),
        lse_ptr
    )
    
    rows = Ref{Cuint}(0)
    cols = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_m_size),
        Cvoid,
        (Ptr{LSE}, Ptr{Cuint}, Ptr{Cuint}),
        lse_ptr, rows, cols
    )
    
    # Convert to Julia array
    return transpose(unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false))
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
    return transpose(unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false))
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
    return transpose(unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false))
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
    return transpose(unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false))
end

function lse_get_ivg_data(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_ivg_data),
        Ptr{ComplexF64},
        (Ptr{LSE},),
        lse_ptr
    )
    
    rows = Ref{Cuint}(0)
    cols = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_ivg_size),
        Cvoid,
        (Ptr{LSE}, Ptr{Cuint}, Ptr{Cuint}),
        lse_ptr, rows, cols
    )
    
    # Convert to Julia array
    return transpose(unsafe_wrap(Array, data_ptr, (rows[], cols[]), own=false))
end

function lse_detImVG(lse_ptr::Ptr{LSE}, E::ComplexF64)
    ccall(Libdl.dlsym(liblse, :lse_detImVG), ComplexF64, (Ptr{LSE}, ComplexF64), lse_ptr, E)
end

function lse_get_psi(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_psi),
        Ptr{ComplexF64},
        (Ptr{LSE},),
        lse_ptr
    )

    rows = Ref{Cuint}(0)
    cols = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_psi_size),
        Cvoid,
        (Ptr{LSE}, Ptr{Cuint}, Ptr{Cuint}),
        lse_ptr, rows, cols
    )
    # Convert to Julia array
    psi = unsafe_wrap(Array, data_ptr, (cols[], rows[]), own=false)
    return transpose(psi)
end

function lse_get_E(lse_ptr::Ptr{LSE})
    data_ptr = ccall(
        Libdl.dlsym(liblse, :lse_get_E),
        Ptr{Cdouble},
        (Ptr{LSE},),
        lse_ptr
    )

    levels = Ref{Cuint}(0)
    
    ccall(
        Libdl.dlsym(liblse, :lse_get_E_size),
        Cvoid,
        (Ptr{Cuint},),
        levels
    )
    # Convert to Julia array
    E = unsafe_wrap(Array, data_ptr, levels[], own=false)
    return E
end

function lse_V(α::Int, β::Int, E::Real, p::ComplexF64, pprime::ComplexF64)
    return ccall(
        Libdl.dlsym(liblse, :V),
        ComplexF64,
        (Cuint, Cuint, Cdouble, ComplexF64, ComplexF64),
        α, β, E, p, pprime
    )
end

function O_00(E::ComplexF64, p::ComplexF64, pprime::ComplexF64, m::Float64)
    return ccall(
        Libdl.dlsym(liblse, :O00),
        ComplexF64,
        (ComplexF64, ComplexF64, ComplexF64, Cdouble),
        E, p, pprime, m
    )
end

function O_01(E::ComplexF64, p::ComplexF64, pprime::ComplexF64, m::Float64)
    return ccall(
        Libdl.dlsym(liblse, :O01),
        ComplexF64,
        (ComplexF64, ComplexF64, ComplexF64, Cdouble),
        E, p, pprime, m
    )
end

function O_10(E::ComplexF64, p::ComplexF64, pprime::ComplexF64, m::Float64)
    return ccall(
        Libdl.dlsym(liblse, :O10),
        ComplexF64,
        (ComplexF64, ComplexF64, ComplexF64, Cdouble),
        E, p, pprime, m
    )
end

function O_11(E::ComplexF64, p::ComplexF64, pprime::ComplexF64, m::Float64)
    return ccall(
        Libdl.dlsym(liblse, :O11),
        ComplexF64,
        (ComplexF64, ComplexF64, ComplexF64, Cdouble),
        E, p, pprime, m
    )
end

export LSE, lse_free, lse_malloc, lse_compute, lse_get_t_data, lse_get_g_data, lse_get_v_data, lse_refresh, lse_gmat, lse_vmat, lse_tmat, lse_get_psi, lse_get_ivg_data, lse_V, lse_detImVG, lse_get_E, O_00, O_01, O_10, O_11, lse_get_iivg_data

end # module
const libwavefunction = joinpath(@__DIR__, "build/linux/x86_64/release/libwavefunction.so")

# Define the function signatures
function wf_new(l::Int, lambda::Real, ngauss::Int)
    ccall(
        (:WFnew, libwavefunction),
        Ptr{Cvoid},
        (UInt64, Float64, UInt64),
        l, lambda, ngauss
    )
end

function wf_free(wf_ptr::Ptr{Cvoid})
    ccall(
        (:WFfree, libwavefunction),
        Cvoid,
        (Ptr{Cvoid},),
        wf_ptr
    )
end

# Get c_solution matrix dimensions
function wf_get_c_solution_dims(wf_ptr::Ptr{Cvoid})
    rows = Ref{Csize_t}(0)
    cols = Ref{Csize_t}(0)
    ccall(
        (:WF_get_c_solution_dims, libwavefunction),
        Cvoid,
        (Ptr{Cvoid}, Ref{Csize_t}, Ref{Csize_t}),
        wf_ptr, rows, cols
    )
    return (rows[], cols[])
end

# Get E_solution vector length
function wf_get_e_solution_length(wf_ptr::Ptr{Cvoid})
    ccall(
        (:WF_get_E_solution_length, libwavefunction),
        Csize_t,
        (Ptr{Cvoid},),
        wf_ptr
    )
end
function wf_get_c_solution_data(wf_ptr::Ptr{Cvoid})
    ccall(
        (:WF_get_c_solution_data, libwavefunction),
        Ptr{Cdouble},
        (Ptr{Cvoid},),
        wf_ptr
    )
end

# Get pointer to E_solution data
function wf_get_e_solution_data(wf_ptr::Ptr{Cvoid})
    ccall(
        (:WF_get_E_solution_data, libwavefunction),
        Ptr{Cdouble},
        (Ptr{Cvoid},),
        wf_ptr
    )
end

# Get tda of c_solution matrix
function wf_get_c_solution_tda(wf_ptr::Ptr{Cvoid})
    ccall(
        (:WF_get_c_solution_tda, libwavefunction),
        Csize_t,
        (Ptr{Cvoid},),
        wf_ptr
    )
end

# Efficient function to get the entire c_solution as a Julia matrix
function get_c_solution_matrix(wf_ptr::Ptr{Cvoid})
    rows, cols = wf_get_c_solution_dims(wf_ptr)
    tda = wf_get_c_solution_tda(wf_ptr)
    data_ptr = wf_get_c_solution_data(wf_ptr)

    # Create a Julia array that references the C memory
    data_array = unsafe_wrap(Array, data_ptr, (tda, cols), own=false)

    # Extract the relevant part (in case tda > rows) and make a copy
    return copy(transpose(data_array[1:rows, 1:cols]))
end

# Efficient function to get the entire E_solution as a Julia vector
function get_e_solution_vector(wf_ptr::Ptr{Cvoid})
    len = wf_get_e_solution_length(wf_ptr)
    data_ptr = wf_get_e_solution_data(wf_ptr)

    # Create a Julia array that references the C memory
    return copy(unsafe_wrap(Array, data_ptr, len, own=false))
end

# Wrapper for psi_n function (real-valued wavefunction in position space)
function psi_n(wf_ptr::Ptr{Cvoid}, r::Float64, n::UInt64, theta::Float64)
    ccall(
        (:psi_n, libwavefunction),
        ComplexF64,
        (Ptr{Cvoid}, Cdouble, UInt64, Cdouble),
        wf_ptr, r, n, theta
    )
end

# Wrapper for psi_n_ft function (complex-valued wavefunction in momentum space)
function psi_n_ft(wf_ptr::Ptr{Cvoid}, p::Float64, n::UInt64)
    ccall(
        (:psi_n_ft, libwavefunction),
        ComplexF64,
        (Ptr{Cvoid}, Cdouble, UInt64),
        wf_ptr, p, n
    )
end

function psi_n_batch(wf_ptr::Ptr{Cvoid}, r_values::Vector{Float64}, n::UInt64, theta::Float64)
    num_points = length(r_values)
    results = Vector{ComplexF64}(undef, num_points)

    ccall(
        (:psi_n_batch, libwavefunction),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{ComplexF64}, Csize_t, UInt64, Cdouble),
        wf_ptr, r_values, results, num_points, n, theta
    )

    return results
end

function psi_n_ft_batch(wf_ptr::Ptr{Cvoid}, p_values::Vector{Float64}, n::UInt64)
    num_points = length(p_values)
    results = Vector{ComplexF64}(undef, num_points)

    ccall(
        (:psi_n_ft_batch, libwavefunction),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{ComplexF64}, Csize_t, UInt64),
        wf_ptr, p_values, results, num_points, n
    )

    return results
end

@inline function quadgauss(f, x::T, w::T) where {T<:Vector{Float64}}
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    @inbounds @simd for i in eachindex(x)
        res += f(x[i]) * w[i]
    end
    return res
end

