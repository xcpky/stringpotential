module LSE

using LinearAlgebra

# Load the shared library
const liblse = joinpath(@__DIR__, "liblse.so")

# Define the LSE struct type (opaque pointer)
mutable struct LSEHandle
    ptr::Ptr{Cvoid}

    # Constructor that ensures proper cleanup
    function LSEHandle(ptr::Ptr{Cvoid})
        handle = new(ptr)
        finalizer(handle) do x
            if x.ptr != C_NULL
                lse_deinit(x)
                x.ptr = C_NULL
            end
        end
        return handle
    end
end

# Function to create a new LSE instance
function new(Ngauss::Integer, Lambda::Real, epsilon::Real)
    ptr = ccall(
        (:lse_new, liblse),
        Ptr{Cvoid},
        (Csize_t, Cdouble, Cdouble),
        Ngauss, Lambda, epsilon
    )

    if ptr == C_NULL
        error("Failed to create LSE instance")
    end

    return LSEHandle(ptr)
end

# Function to run the LSE solver
function compute(lse::LSEHandle, E::Real)
    result = ccall(
        (:lse_compute, liblse),
        Cint,
        (Ptr{Cvoid}, Cdouble),
        lse.ptr, E
    )

    if result != 0
        error("LSE solver failed with error code $result")
    end

    return nothing
end

# Function to clean up resources
function lse_deinit(lse::LSEHandle)
    if lse.ptr != C_NULL
        ccall(
            (:lse_deinit, liblse),
            Cvoid,
            (Ptr{Cvoid},),
            lse.ptr
        )
    end
    return nothing
end

# Helper function to get matrix dimensions
function get_g_matrix_size(lse::LSEHandle)
    rows_ref = Ref{Cuint}(0)
    cols_ref = Ref{Cuint}(0)

    ccall(
        (:lse_get_g_size, liblse),
        Cvoid,
        (Ptr{Cvoid}, Ref{Cuint}, Ref{Cuint}),
        lse.ptr, rows_ref, cols_ref
    )

    return (rows_ref[], cols_ref[])
end

function get_v_matrix_size(lse::LSEHandle)
    rows_ref = Ref{Cuint}(0)
    cols_ref = Ref{Cuint}(0)

    ccall(
        (:lse_get_v_size, liblse),
        Cvoid,
        (Ptr{Cvoid}, Ref{Cuint}, Ref{Cuint}),
        lse.ptr, rows_ref, cols_ref
    )

    return (rows_ref[], cols_ref[])
end

function get_t_matrix_size(lse::LSEHandle)
    rows_ref = Ref{Cuint}(0)
    cols_ref = Ref{Cuint}(0)

    ccall(
        (:lse_get_t_size, liblse),
        Cvoid,
        (Ptr{Cvoid}, Ref{Cuint}, Ref{Cuint}),
        lse.ptr, rows_ref, cols_ref
    )

    return (rows_ref[], cols_ref[])
end

function gmat(lse::LSEHandle, E::Real)
    result = ccall(
        (:lse_gmat, liblse),
        Cint,
        (Ptr{Cvoid}, Cdouble),
        lse.ptr, E
    )
    if result != 0
        error("LSE G matrix failed with error code $result")
    end
    return get_g_matrix(lse)
end

# Function to get G matrix
function get_g_matrix(lse::LSEHandle)
    rows, cols = get_g_matrix_size(lse)

    # Get pointer to the raw data
    data_ptr = ccall(
        (:lse_get_g_data, liblse),
        Ptr{Cdouble},
        (Ptr{Cvoid},),
        lse.ptr
    )

    # Complex data is stored as alternating real and imaginary parts
    # So the actual array size is 2 * rows * cols
    data = unsafe_wrap(Array, data_ptr, 2*rows*cols, own=false)

    # Create a complex matrix from the raw data
    return reshape(ComplexF64.(data[1:2:end], data[2:2:end]), (rows, cols))
end

# Function to get V matrix
function get_v_matrix(lse::LSEHandle)
    rows, cols = get_v_matrix_size(lse)

    # Get pointer to the raw data
    data_ptr = ccall(
        (:lse_get_v_data, liblse),
        Ptr{Cdouble},
        (Ptr{Cvoid},),
        lse.ptr
    )

    # Complex data is stored as alternating real and imaginary parts
    data = unsafe_wrap(Array, data_ptr, 2*rows*cols, own=false)

    # Create a complex matrix from the raw data
    return reshape(ComplexF64.(data[1:2:end], data[2:2:end]), (rows, cols))
end

# Function to get T matrix
function get_t_matrix(lse::LSEHandle)
    rows, cols = get_t_matrix_size(lse)

    # Get pointer to the raw data
    data_ptr = ccall(
        (:lse_get_t_data, liblse),
        Ptr{Cdouble},
        (Ptr{Cvoid},),
        lse.ptr
    )

    # Complex data is stored as alternating real and imaginary parts
    data = unsafe_wrap(Array, data_ptr, 2*rows*cols, own=false)

    # Create a complex matrix from the raw data
    return reshape(ComplexF64.(data[1:2:end], data[2:2:end]), (rows, cols))
end

# Export the main functions
export new, run, compute, get_g_matrix, get_v_matrix, get_t_matrix, get_g_matrix, gmat

end # module
