include("lse.jl")

Ngauss = 48
E = delta[1]-2:0.002001:delta[2]+0.5

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
    tmatrices = tmat.(4, E, Ngauss)
    osT = onshellT.(tmatrices)
    len = size(E)[1]
    T = [[osT[i][1] for i in 1:len], [osT[i][2] for i in 1:len], [osT[i][1] for i in 1:len], [osT[i][4] for i in 1:len]]
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
end
