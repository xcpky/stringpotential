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

@inline function Delta011(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_11(p, pprime) - E
    D = ωprime_11(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return 1 / B * log(Complex((b + C) * (b + D) / (a + C) / (a + D)))
end

@inline function Delta012(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_12(p, pprime) - E
    D = ωprime_12(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return 1 / B * log(Complex((b + C) * (b + D) / (a + C) / (a + D)))
end

@inline function Delta021(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_21(p, pprime) - E
    D = ωprime_21(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return 1 / B * log(Complex((b + C) * (b + D) / (a + C) / (a + D)))
end

@inline function Delta022(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_11(p, pprime) - E
    D = ωprime_11(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return 1 / B * log(Complex((b + C) * (b + D) / (a + C) / (a + D)))
end

@inline function Delta111(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_11(p, pprime) - E
    D = ωprime_11(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return -((C + D) * (a - b) + 2 * B + (A - C * C) * log(Complex((a + C) / (b + C))) + (A - D * D) * log(Complex((a + D) / (b + D)))) / B / B
end

@inline function Delta112(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_12(p, pprime) - E
    D = ωprime_12(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return -((C + D) * (a - b) + 2 * B + (A - C * C) * log(Complex((a + C) / (b + C))) + (A - D * D) * log(Complex((a + D) / (b + D)))) / B / B
end

@inline function Delta121(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_21(p, pprime) - E
    D = ωprime_21(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return -((C + D) * (a - b) + 2 * B + (A - C * C) * log(Complex((a + C) / (b + C))) + (A - D * D) * log(Complex((a + D) / (b + D)))) / B / B
end

@inline function Delta122(E, p, pprime, m)
    A = p^2 + pprime^2 + m^2
    B = 2 * p * pprime
    C = ω_22(p, pprime) - E
    D = ωprime_22(p, pprime) - E
    a = sqrt(Complex(A - B))
    b = sqrt(Complex(A + B))
    return -((C + D) * (a - b) + 2 * B + (A - C * C) * log(Complex((a + C) / (b + C))) + (A - D * D) * log(Complex((a + D) / (b + D)))) / B / B
end

@inline function O11(E, p, pprime, m)
    -3 * g_pi * g_pi / 24 / f_pi / f_pi * (2 * p * pprime * Delta111(E, p, pprime, m) - (p^2 + pprime^2) * Delta011(E, p, pprime, m))
end

@inline function O12(E, p, pprime, m)
    -3 * g_pi * g_pi / 24 / f_pi / f_pi * (2 * p * pprime * Delta112(E, p, pprime, m) - (p^2 + pprime^2) * Delta012(E, p, pprime, m))
end

@inline function O21(E, p, pprime, m)
    -3 * g_pi * g_pi / 24 / f_pi / f_pi * (2 * p * pprime * Delta121(E, p, pprime, m) - (p^2 + pprime^2) * Delta021(E, p, pprime, m))
end

@inline function O22(E, p, pprime, m)
    -3 * g_pi * g_pi / 24 / f_pi / f_pi * (2 * p * pprime * Delta122(E, p, pprime, m) - (p^2 + pprime^2) * Delta022(E, p, pprime, m))
end

@inline function V_OME_11(E, p, pprime)
    return 3 * O11(E, p, pprime, m_pi) + O11(E, p, pprime, m_eta)
end

@inline function V_OME_12(E, p, pprime)
    return 2 * sqrt(2) * O12(E, p, pprime, m_K)
end

@inline function V_OME_21(E, p, pprime)
    return 2 * sqrt(2) * O21(E, p, pprime, m_K)
end

@inline function V_OME_22(E, p, pprime)
    return 4 / 3 * O22(E, p, pprime, m_eta)
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

function V(α, β, E, p, pprime)
    if α == 0 && β == 0
        V_OME_11(E, p, pprime) 
    elseif α == 0 && β == 1
        V_OME_12(E, p, pprime) 
    elseif α == 1 && β == 0
        V_OME_21(E, p, pprime)
    else
        V_OME_22(E, p, pprime)
    end
end


