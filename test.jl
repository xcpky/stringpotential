using CSV; using DataFrames; ENV["COLUMNS"]=160
using Plots
using Interpolations
using QuadGK; using LaTeXStrings
# using BenchmarkTools
using LinearAlgebra, StaticArrays
include("hadronmasses.jl")
include("quadgauss.jl")
const C1=0.2005;
const C2=-0.0657;
const ϵ=1e-10;
const XW = Ref(gauss(64, 2mπ+0.000001, 5));
const xx1, ww1 = gauss(32, -1, 1);
const xx2, ww2 = gauss(32, 0, 1);
σ(s, m) = sqrt(1 - 4m^2/s)
λ(x, y, z) = x^2 + y^2 + z^2 -2x*y - 2y*z - 2z*x
q3(s) = sqrt(λ(mψp^2, mjψ^2, s))/(2mψp)

function qsq(E,m1,m2) 
    return λ(E^2,m1^2,m2^2)/(4E^2)
end
function xsqrt(x)
    imag(x) >=0 ? sqrt(x+0im) : -sqrt(x-0im)
end

# here only consider the S-wave part, par = [c1, c2]
function amp0(s, par;which=:BO)
    c1, c2 = par
    msq = mψp^2; mpisq = mπ^2
    qsq = λ(msq, mjψ^2, s)/(4msq)
    σsq = 1 - 4mpisq/s
    return -2/fπ^2 * (c1 * (s-2mpisq) 
        + 0.5c2 * (s + qsq * (1 - σsq/3)) ) * absomnes(sqrt(s),which=which)
end
function vσGau1s(p1,p2, par, Λ = 1.0;which=:BO)
    function vintegrand(w)
        s = w^2
        sqrt(s-4mπ^2) * (amp0(s, par,which=which))^2*exp(-s/Λ^2)*log((s+(p1+p2)^2)/(s+(p1-p2)^2))
    end
    xxv, wwv=XW[]
    res = -1/(128π^2*mjψ^2*p1*p2) *  quadgauss(vintegrand, xxv, wwv) 
    return res
end
function vσHC(p1,p2, par, Λ = 1.0;which=:BO)
    function vintegrand(w)
        s = w^2
        sqrt(s-4mπ^2) * (amp0(s, par,which=which))^2*log((s+(p1+p2)^2)/(s+(p1-p2)^2))
    end
     vintegrand1(x)=(Λ-2mπ)*vintegrand((Λ-2mπ)*x+2mπ)

     res = -1/(128π^2*mjψ^2*p1*p2) *  quadgauss(vintegrand1, xx2, ww2) 
    return res
end
function Vmat(E,Λ,Vcont=0,N=16;which=:BO)
    pole = 6.0
    if abs(E-pole) < 1e-8
        E += 1e-6
    end
    xx, ww = gauss(N,0, Λ)
    xx=convert(Array{Complex},xx)
    ww=convert(Array{Complex},ww)
    x0=xsqrt(qsq(E,mjψ,mjψ))
    append!(xx,x0)
    Vmatrix = zeros(ComplexF64,N+1,N+1);
    for i=1:1:N+1
        for j=1:1:N+1
            Vmatrix[i,j]=0.001/(E-pole)
        end
    end
    return Vmatrix
end
function Gmat(E,Λ,rs=1,N=16)
    dE=E-2mjψ
    μ=mjψ/2
    x0=0
    if rs==1
        x0=xsqrt(2μ*dE)
    elseif rs==2
        x0=-xsqrt(2μ*dE)
    end
    Gmatrix = zeros(ComplexF64,N+1,N+1);
    xx, ww = gauss(N,0, Λ)
    integrand(q)=1/(dE-q^2/2μ+im*ϵ)
    Gmatrix[N+1,N+1]=1/(2*π^2)*x0*(μ*(-im*π + log((Λ+x0)/(Λ-x0)))
        -x0*quadgauss(integrand,xx,ww))
    for i=1:N
        Gmatrix[i,i]=1/(2*π^2)*xx[i]^2*ww[i]/(dE-xx[i]^2/2μ+im*ϵ)
    end
    return Gmatrix
end
function T(E,Λ,Vcont=0,N=16;which=:BO)
    V=Vmat(E,Λ,Vcont,N;which)
    G=Gmat(E,Λ,1,N)
    Tv=inv(I-V*G)*V[:,N+1]
    return Tv[end]
end

E = 5.5:0.01:6.8
ot = T.(E, 2, 0, 64)
plot(E, abs.(ot),dpi=300)
savefig("test.png")
