include("constants.jl")
@inline function xsqrt(x)
    imag(x) >= 0 ? sqrt(x + 0im) : -sqrt(x - 0im)
end

@inline function ω_11(p, pprime)
    2 * m_B + (p^2 + pprime^2) / 2 / m_B
end

@inline function ω_12(p, pprime)
    m_B + pprime^2 / 2 / m_B + m_B_s + p^2 / 2 / m_B_s
end

@inline function ω_21(p, pprime)
    m_B_s + pprime^2 / 2 / m_B_s + m_B + p^2 / 2 / m_B
end

@inline function ω_22(p, pprime)
    2 * m_B_s + (p^2 + pprime^2) / 2 / m_B_s
end

@inline function ωprime_11(p, pprime)
    2 * m_B_star + (p^2 + pprime^2) / 2 / m_B_star
end

@inline function ωprime_12(p, pprime)
    m_B_star + pprime^2 / 2 / m_B_star + m_B_star_s + p^2 / 2 / m_B_star_s
end

@inline function ωprime_21(p, pprime)
    m_B_star_s + pprime^2 / 2 / m_B_star_s + m_B_star + p^2 / 2 / m_B_star
end

@inline function ωprime_22(p, pprime)
    2 * m_B_star_s + (p^2 + pprime^2) / 2 / m_B_star_s
end

const ω = (
    (ω_11, ω_12),
    (ω_21, ω_22)
)

const ωprime = (
    (ωprime_11, ωprime_12),
    (ωprime_21, ωprime_22)
)

@inline function O(::Val{α}, ::Val{β}, E, p, pprime, m) where {α,β}
    return -1 / 4 / p / pprime * (log(Complex((E - (m + (p - pprime)^2 / 2 / m) - ω[α][β](p, pprime)) / E - (m + (p + pprime)^2 / 2 / m) - ω[α][β](p, pprime))) + log(Complex((E - (m + (p - pprime)^2 / 2 / m) - ωprime[α][β](p, pprime)) / E - (m + (p + pprime)^2 / 2 / m) - ωprime[α][β](p, pprime))))
end

@inline function V_OME_11(E, p, pprime)
    -3 * (3 * O(Val(1), Val(1), E, p, pprime, m_pi) + O(Val(1), Val(1), E, p, pprime, m_eta) / 3)
end

@inline function V_OME_12(E, p, pprime)
    2^(3 / 2) * O(Val(1), Val(2), E, p, pprime, m_K)
end

@inline function V_OME_21(E, p, pprime)
    2^(3 / 2) * O(Val(2), Val(1), E, p, pprime, m_K)
end

@inline function V_OME_22(E, p, pprime)
    2 / 3 * O(Val(2), Val(2), E, p, pprime, m_eta)
end

@inline function Ctct_11(g_C)
    -3 * 2 * g_C
end

@inline function Ctct_12(g_C)
    (1 / 2)^(3 / 2) * 4 * g_C
end

@inline function Ctct_22(g_C)
    g_C
end
