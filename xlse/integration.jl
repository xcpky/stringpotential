using QuadGK
xsqrtright(x) = imag(x) >= 0 ? sqrt(x + 0im) : -sqrt(x - 0im)
xsqrtup(x) = imag(x) >= 0 && real(x) < 0 ? -sqrt(x + 0im) : sqrt(x - 0im)
xsqrtleft(x) = sqrt(Complex(x + 0im))
xsqrtdn(x) = imag(x) < 0 && real(x) < 0 ? -sqrt(x - 0im) : sqrt(x + 0im)

xsqrt = xsqrtdn
epsi = -1e-6
B = im

f1(x) = 1 / xsqrt(epsi - x * im) / (1 + xsqrt(epsi - x * im))
f2(x) = 1 / xsqrt(-epsi - x * im) / (1 + xsqrt(-epsi - x * im))
f(A, B, x) = 1 / xsqrt(A - B * x) / (1 + xsqrt(A - B * x))
inte(a) = quadgk(x -> f(a, B, x), -1, 1)[1]
I(A, B) = 2 / B * (log(xsqrt(A + B) + 1) - log(xsqrt(A - B) + 1))
I(A, B) = 2 / B * (log((xsqrt(A + B) + 1) / (xsqrt(A - B) + 1)))

F1(x) = -2 / im * log(xsqrt(1 - x * im) + 1)
F2(x) = -2 / im * log(xsqrt(-1 - x * im) + 1)

x, w = gauss(128, 0, 1)

path = -1 .+ x .* 2;
weight = w .* 2;

pathn = -1 .+ x .* 1;
weightn = w .* 1;

pathp = 0 .+ x .* 1;
weightp = weightn

pathsegs = [pathn; pathp]
weightsegs = [weightn; weightp]

a = -1.2im
pathl = [-1 .+ x .* a; (-1 + a) .+ x .* 2; (1 + a) .+ x .* -a]
weightl = [w .* a; w .* 2; w .* -a]

a = 1.2im
pathr = [-1 .+ x .* a; (-1 + a) .+ x .* 2; (1 + a) .+ x .* -a]
weightr = [w .* a; w .* 2; w .* -a]

val1 = []
push!(val1, quadgk(f1, -1, 1)[1])
tmp = 0
for i in eachindex(path)
    global tmp += f1(path[i]) * weight[i]
end
push!(val1, tmp)

push!(val1, F1(1) - F1(-1))

val2 = []
push!(val2, quadgk(f2, -1, 1)[1])
tmp = 0;
for i in eachindex(path)
    global tmp += f2(path[i]) * weight[i]
end
push!(val2, tmp)

tmp = 0
for i in eachindex(pathsegs)
    global tmp += f2(pathsegs[i]) * weightsegs[i]
end
push!(val2, tmp)

tmp = 0
for i in eachindex(pathr)
    global tmp += f2(pathr[i]) * weightr[i]
end
push!(val2, tmp)

push!(val2, F2(1) - F2(-1))


xsqrt = xsqrtright


tmp = 0
for i in eachindex(pathsegs)
    global tmp += f1(pathsegs[i]) * weightsegs[i]
end
push!(val1, tmp)

tmp = 0
for i in eachindex(pathl)
    global tmp += f1(pathl[i]) * weightl[i]
end
push!(val1, tmp)

push!(val1, F1(1) - F1(-1))


push!(val2, quadgk(f2, -1, 1)[1])

tmp = 0
for i in eachindex(path)
    global tmp += f2(path[i]) * weight[i]
end

push!(val2, tmp)

push!(val2, F2(1) - F2(-1))

xsqrt = xsqrtdn


push!(val1, quadgk(f1, -1, 1)[1])

tmp = 0
for i in eachindex(path)
    global tmp += f1(path[i]) * weight[i]
end
push!(val1, tmp)

push!(val1, F1(1) - F1(-1))



push!(val2, quadgk(f2, -1, 1)[1])

tmp = 0
for i in eachindex(path)
    global tmp += f2(path[i]) * weight[i]
end
push!(val2, tmp)

push!(val2, F2(1) - F2(-1))
