include("lse.jl")

Lambda = 4
pNgauss = 64
# E = delta[1]-2:0.002001:delta[2]+0.5
Erange = LinRange(-0.9, delta[1] + 0.4, 500)
Erange = LinRange(-0.5, 0.5, 500)
# E = -1.5:0.00201:0.4

function onshellG(matrix)
    return [tr(matrix[1:pNgauss+1, 1:pNgauss+1]), tr(matrix[pNgauss+2:2*pNgauss+2, pNgauss+2:2*pNgauss+2])]
end

function onshellT(matrix)
    return [matrix[pNgauss+1, pNgauss+1], matrix[pNgauss+1, end], matrix[end, pNgauss+1], matrix[end, end]]
end

if "--onshellG" in ARGS
    gmatrices = gmat.(4, Erange, pNgauss)
    tce = onshellG.(gmatrices)
    g11 = [tce[i][1] for i in 1:size(Erange)[1]]
    g22 = [tce[i][2] for i in 1:size(Erange)[1]]
    using Plots
    plot(Erange, real.(g11), dpi=300)
    plot!(Erange, imag.(g11))
    plot!(Erange, real.(g22))
    plot!(Erange, imag.(g22))
    ylims!(-1.5, 0.5)
    savefig("onshellG.png")
end

if "--onshellT" in ARGS
    osT = onshellT.(tmat.(Lambda, Erange, pNgauss))
    len = size(Erange)[1]
    T = [[abs(osT[i][1]) for i in 1:len], [abs(osT[i][2]) for i in 1:len], [abs(osT[i][1]) for i in 1:len], [abs(osT[i][4]) for i in 1:len]]
    using Plots
    plot(Erange, T[1], dpi=300)
    for i in 2:4
        plot!(Erange, T[i])
    end
    vline!(delta, ls=:dash)
    # ylims!(0, 10)
    # ylims!(0, 1e5)
    savefig("onshellT.png")
end

if "--onshellV" in ARGS
    osT = onshellT.(vmat.(Lambda, Erange, pNgauss, 0))
    len = size(Erange)[1]
    T = [[abs(osT[i][1]) for i in 1:len], [abs(osT[i][2]) for i in 1:len], [abs(osT[i][3]) for i in 1:len], [abs(osT[i][4]) for i in 1:len]]
    using Plots
	plot(dpi=400)
    plot!(Erange, T[1])
	for i in [2, 3, 4]
        plot!(Erange, T[i])
    end
    vline!(delta, ls=:dash)
    # ylims!(0, 10)
    # ylims!(0, 1e5)
    savefig("onshellV.png")
end

if "--Det" in ARGS
    de = detImVG.(Lambda, Erange, pNgauss)
    len = size(de)[1]
    using Plots
    plot(Erange, abs.(de), dpi=400)
    # ylims!(0, 100)
    savefig("det.png")

end

if "--test" in ARGS
    pNgauss = 3
    Lambda = 4
    epsilon = 1e-9
    Erange = -0.3
    n = 2 * (pNgauss + 1)

    using Libdl

    # Load the shared library
    const libscript = Libdl.dlopen(joinpath(@__DIR__, "xlse/build/linux/x86_64/release/libscript.so"))
    lse = ccall(Libdl.dlsym(libscript, :lse_malloc), Ptr{Cvoid}, (Csize_t, Cdouble, Cdouble), pNgauss, Lambda, epsilon)
    ccall(Libdl.dlsym(libscript, :lse_refresh), Cvoid, (Ptr{Cvoid}, ComplexF64, Ptr{Cdouble}, Cuint), lse, Erange, [0.0, 0, 0, 0], 3)
    ccall(Libdl.dlsym(libscript, :lse_vmat), Cint, (Ptr{Cvoid},), lse)
    ptr = ccall(Libdl.dlsym(libscript, :lse_get_v_data), Ptr{ComplexF64}, (Ptr{Cvoid},), lse)
    vv = transpose(copy(unsafe_wrap(Array, ptr, (n, n), own=false)))


    p = xsqrt(2 * mu[1] * (Erange - delta[1]))
    data = vmat(Lambda, Erange, pNgauss)
    Erange += m11 + m12
    function testfunc(E, p, pprime, m)
        log(Complex((E - (m + (p - pprime)^2 / 2 / m) - ω[1][1](p, pprime)) / (E - (m + (p + pprime)^2 / 2 / m) - ω[1][1](p, pprime))))
    end
    xi, wi = gauss(pNgauss, 0, Lambda)
    m = m_pi
    pprime = p
    a1 = (Erange - (m + (p - pprime)^2 / 2 / m) - ω[1][1](p, pprime))
    b1 = (Erange - (m + (p + pprime)^2 / 2 / m) - ω[1][1](p, pprime))
    a2 = Erange - (m + (p - pprime)^2 / 2 / m) - ωprime[1][1](p, pprime)
    b2 = Erange - (m + (p + pprime)^2 / 2 / m) - ωprime[1][1](p, pprime)
    part1 = log(Complex(a1 / b1))
    part2 = log(Complex(a2 / b2))
    ooo = -g_pi^2 / f_pi^2 / p / pprime * (part1 + part2)


end

if "--testSchrodinger" in ARGS
    include("schrodingerEquation.jl")
    using .SchrodingerEquation
    Ngauss = 800
    Nlag = 5
    Scaling = 1
    p, w = gauss(Ngauss)
    @. p = Scaling * ((1 + p) / (1 - p))
    @. w = Scaling * (2 * w / (1 - p)^2)
    y = @. (p^2 + p'^2) / 2 / (p * p')
    sol = Solver(1.0, 1.0, 1.0, 0.0)
    Erange, psi = solve(sol, 0)
end

if "--integration" in ARGS
    Erange = LinRange(0.1, 0.8, 256)
    println(pNgauss)
    integrand(x, E, p1) = V_OME_11(E, p1, x) * x^2 / 2 / pi / 2 / (E - m11 - m12 - x^2 / 2 / mu[1] + ϵ * im)
    p, w = gauss(pNgauss, 0, Lambda)
    res = Array{ComplexF64}(undef, pNgauss, length(Erange))
    for i in 1:pNgauss
        for j in eachindex(Erange)
            res[i, j], _ = quadgk(x -> integrand(x, Erange[j] + m11 + m12, p[i]), 0, Lambda)
        end
    end
    using Plots
    plot(dpi=400)
    surface!(Erange, p, abs.(res))
    savefig("quadgk.png")
end

function naivedif1(Ngauss, Nlag, p)
    dif = zeros(Float64, (Ngauss, Ngauss))
    for j in 1:Ngauss
        for i in 1:Ngauss
            start = min(Ngauss - Nlag + 1, max(1, i - Int(floor(Nlag // 2))))
            range = start:start+Nlag-1
            mask = setdiff(range, j)
            for m in mask
                mtemp = 1
                mas = setdiff(mask, m)
                for k in mas
                    mtemp *= (p[i] - p[k]) / (p[j] - p[k])
                end
                dif[i, j] += mtemp / (p[j] - p[m])
            end
        end
    end
    return dif
end

function naivedif2(Ngauss, Nlag, p)
    dif = zeros(Float64, (Ngauss, Ngauss))
    for j in 1:Ngauss
        for i in 1:Ngauss
            start = min(Ngauss - Nlag + 1, max(1, i - Int(floor(Nlag // 2))))
            range = start:start+Nlag-1
            mask = setdiff(range, j)
            for m in mask
                mas = setdiff(mask, m)
                mtemp = 0
                for l in mas
                    ma = setdiff(mas, l)
                    ltemp = 1
                    for k in ma
                        ltemp *= (p[i] - p[k]) / (p[j] - p[k])
                    end
                    mtemp += ltemp / (p[j] - p[l])
                end
                dif[i, j] += mtemp / (p[j] - p[m])
            end
        end
    end
    return dif
end

function naivedif(Ngauss, Nlag, p)
    return naivedif1(Ngauss, Nlag, p), naivedif2(Ngauss, Nlag, p)
end
if "--REPL" in ARGS
    local E = -0.1
    p = 0.2
    while true
        pprime = parse(Float64, readline())
        if abs(pprime) < 1e-8
            exit(0)
        end
        # v = Vπu(E, p, pprime, m_B, 0, m_B_star, gamma_B_star, m_B, 0, m_B_star, gamma_B_star, m_pi, 1)
        v = V22(E + m11 + m12, p, pprime )
        # v = onshellT(vmat(Lambda, E, pNgauss, 0))[1]
        println(v)
    end
end

