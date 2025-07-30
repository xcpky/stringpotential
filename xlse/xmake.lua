add_rules("mode.debug", "mode.release")
set_languages("gnu23")
add_cflags("-Wall", "-Wextra", "-DHAVE_INLINE", "-fgnuc-version=8", "-Wconversion")
set_toolchains("clang")
add_requires("gsl")
-- add_requires("matplotplusplus")
add_rules("plugin.compile_commands.autoupdate")

target("wavefunction")
    set_kind("shared")
    add_files("src/wavefunction.c", "src/constants.c")
    add_packages("gsl")
    add_links("m")

target("lse")
    set_kind("shared")
    -- add_headerfiles("src/lse.h")
    add_files("src/lse.c", "src/constants.c", "src/wavefunction.c", "src/ome.c")
    add_packages("gsl")
    add_links("m")

target("xlse")
    set_kind("binary")
    add_deps("lse")
    add_deps("wavefunction")
    add_deps("script")
    add_packages("gsl")
    add_links("m")
    add_files("src/main.c", "src/constants.c")

target("script")
    local cpu_count = os.cpuinfo().ncpu
    add_defines("NTHREADS=" .. cpu_count)
    set_kind("shared")
    add_files("src/script.c", "src/constants.c", "src/wavefunction.c", "src/lse.c", "src/ome.c")
    add_packages("gsl")
    add_links("m")
--
-- target("xlsepp")
--     set_kind("binary")
--     add_deps("lse")
--     set_languages("c++20")
--     add_deps("wavefunction")
--     add_deps("script")
--     add_packages("matplotplusplus")
--     add_packages("gsl")
--     add_links("m")
--     add_files("src/main.cpp", "src/constants.c")
--
--
-- If you want to known more usage about xmake, please see https://xmake.io
--
-- ## FAQ
--
-- You can enter the project directory firstly before building project.
--
--   $ cd projectdir
--
-- 1. How to build project?
--
--   $ xmake
--
-- 2. How to configure project?
--
--   $ xmake f -p [macosx|linux|iphoneos ..] -a [x86_64|i386|arm64 ..] -m [debug|release]
--
-- 3. Where is the build output directory?
--
--   The default output directory is `./build` and you can configure the output directory.
--
--   $ xmake f -o outputdir
--   $ xmake
--
-- 4. How to run and debug target after building project?
--
--   $ xmake run [targetname]
--   $ xmake run -d [targetname]
--
-- 5. How to install target to the system directory or other output directory?
--
--   $ xmake install
--   $ xmake install -o installdir
--
-- 6. Add some frequently-used compilation flags in xmake.lua
--
-- @code
--    -- add debug and release modes
--    add_rules("mode.debug", "mode.release")
--
--    -- add macro definition
--    add_defines("NDEBUG", "_GNU_SOURCE=1")
--
--    -- set warning all as error
--    set_warnings("all", "error")
--
--    -- set language: c99, c++11
--    set_languages("c99", "c++11")
--
--    -- set optimization: none, faster, fastest, smallest
--    set_optimize("fastest")
--
--    -- add include search directories
--    add_includedirs("/usr/include", "/usr/local/include")
--
--    -- add link libraries and search directories
--    add_links("tbox")
--    add_linkdirs("/usr/local/lib", "/usr/lib")
--
--    -- add system link libraries
--    add_syslinks("z", "pthread")
--
--    -- add compilation and link flags
--    add_cxflags("-stdnolib", "-fno-strict-aliasing")
--    add_ldflags("-L/usr/local/lib", "-lpthread", {force = true})
--
-- @endcode
--

