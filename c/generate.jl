using Clang.Generators
using Clang.LibClang.Clang_jll

# Path to your C header files
headers = ["lse.h"]

# Include directories (add all necessary include paths)
include_dirs = [
    "/usr/include/gsl",
    "/usr/include",
    "/usr/local/include",
]


# Create options for the generator
options = load_options(joinpath(@__DIR__, "generator.toml"))

# If you don't have a generator.toml file, you can set options directly:
# options = Dict{String,Any}(
#     "output_file_path" => output_dir,
#     "library_name" => "YourLibrary",
#     "output_mode" => "single",  # or "multiple" for multiple files
#     "header_wrapped" => (header, cursor) -> true,  # Include all headers
#     "header_library" => header -> "YourLibrary",
#     "clang_args" => vcat(map(dir -> "-I$dir", include_dirs)),
#     "clang_diagnostics" => true
# )

# Create the context
ctx = create_context(headers, options["general"]["clang_args"], options)

# Build the bindings
build!(ctx)
