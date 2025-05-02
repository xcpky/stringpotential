# Load the shared library
include("lse.jl")
using .CLSE
using QuadGK
using LinearAlgebra
const libwavefunction = joinpath(@__DIR__, "libwavefunction.so")

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
    return copy(data_array[1:rows, 1:cols])
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

# Example usage
# wf = wf_new(0, 20, 40)
# c_matrix = get_c_solution_matrix(wf)
# e_vector = get_e_solution_vector(wf)
# println("C solution matrix: ", c_matrix)
Ngauss = 100
lse = lse_malloc(Ngauss, 4, 1e-7)
lse_refresh(lse, -0.2)
function OnshellT(lse::Ptr{LSE}, E::Real)
    lse_compute(lse, E)
    T = lse_get_t_data(lse)
    return T[Ngauss+1, Ngauss+1], T[Ngauss+1, end], T[end, Ngauss+1], T[end]
end
function OnshellG(lse::Ptr{LSE}, E::Real)
    lse_refresh(lse, E)
    lse_gmat(lse)
    G = lse_get_g_data(lse)
    return [tr(G[1:Ngauss+1, 1:Ngauss+1]), tr(G[Ngauss+2:2*Ngauss+2, Ngauss+2:2*Ngauss+2])]
end
function OnshellV(lse::Ptr{LSE}, E::Real)
    lse_refresh(lse, E)
    lse_vmat(lse)
    V = lse_get_v_data(lse)
    return V[Ngauss+1, Ngauss+1], V[Ngauss+1, end], V[end, Ngauss+1], V[end]
end
function OnshellPsi(lse::Ptr{LSE}, E::Real, n::Int)
    lse_refresh(lse, E)
    psi = lse_get_psi(lse)
    return psi[n, end-2], psi[n, end-1]
end

E = -1:0.02:1
if "--onshellPsi" in ARGS
    op = OnshellPsi.(lse, E, 1)
    psi0 = [op[i][1] for i in eachindex(E)]
    psi1 = [op[i][2] for i in eachindex(E)]
end

if "--onshellG" in ARGS
    tce = OnshellG.(lse, E)
    g11 = [tce[i][1] for i in eachindex(E)]
    g22 = [tce[i][2] for i in eachindex(E)]
end

if "--onshellV" in ARGS
    oV = OnshellV.(lse, E)
    oV11 = [abs(oV[i][1]) for i in eachindex(E)]
    oV12 = [abs(oV[i][2]) for i in eachindex(E)]
    oV21 = [abs(oV[i][3]) for i in eachindex(E)]
    oV22 = [abs(oV[i][4]) for i in eachindex(E)]
end

if "--onshellT" in ARGS
    oT = OnshellT.(lse, E)
    oT11 = [abs(oT[i][1]) for i in eachindex(E)]
    oT12 = [abs(oT[i][2]) for i in eachindex(E)]
    oT21 = [abs(oT[i][3]) for i in eachindex(E)]
    oT22 = [abs(oT[i][4]) for i in eachindex(E)]
end

