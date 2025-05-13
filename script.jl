include("lse.jl")

Ngauss = 40
# E = delta[1]-2:0.002001:delta[2]+0.5
E = -2.5:0.00201:0.5

function onshellG(matrix)
    return [tr(matrix[1:Ngauss+1, 1:Ngauss+1]), tr(matrix[Ngauss+2:2*Ngauss+2, Ngauss+2:2*Ngauss+2])]
end

function onshellT(matrix)
    return [matrix[Ngauss+1, Ngauss+1], matrix[Ngauss+1, end], matrix[end, Ngauss+1], matrix[end, end]]
end

if "--onshellG" in ARGS
    gmatrices = gmat.(4, E, Ngauss)
    tce = onshellG.(gmatrices)
    g11 = [tce[i][1] for i in 1:size(E)[1]]
    g22 = [tce[i][2] for i in 1:size(E)[1]]
end

if "--onshellT" in ARGS
    osT = onshellT.(gmat.(4, E, Ngauss))
    len = size(E)[1]
    T = [[abs(osT[i][1]) for i in 1:len], [abs(osT[i][2]) for i in 1:len], [abs(osT[i][1]) for i in 1:len], [abs(osT[i][4]) for i in 1:len]]
    using Plots
    plot(E, T[1], dpi=300)
    for i in 2:4
        plot!(E, T[i])
    end
    # ylims!(0, 1e5)
    savefig("onshellT.png")
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
    E, psi = solve(sol, 0)
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
