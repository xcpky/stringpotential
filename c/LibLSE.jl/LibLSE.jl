module cLSE

struct WaveFunction
    l::Cint
    Lambda::Cdouble
    Ngauss::Cint
    xi::Ptr{Cdouble}
    wi::Ptr{Cdouble}
    table::Ptr{Cint}
    c_solution::Ptr{Cint}
    E_solution::Ptr{Cint}
end

function WFnew(l, Lambda, Ngauss)
    @ccall LibLSE.WFnew(l::Cint, Lambda::Cdouble, Ngauss::Cint)::Ptr{WaveFunction}
end

function WFfree(self)
    @ccall LibLSE.WFfree(self::Ptr{WaveFunction})::Cvoid
end

function build(self)
    @ccall LibLSE.build(self::Ptr{WaveFunction})::Cvoid
end

function WF_get_c_solution_dims(self, rows, cols)
    @ccall LibLSE.WF_get_c_solution_dims(self::Ptr{WaveFunction}, rows::Ptr{Cint}, cols::Ptr{Cint})::Cvoid
end

function WF_get_c_solution_data(self)
    @ccall LibLSE.WF_get_c_solution_data(self::Ptr{WaveFunction})::Ptr{Cdouble}
end

function WF_get_c_solution_tda(self)
    @ccall LibLSE.WF_get_c_solution_tda(self::Ptr{WaveFunction})::Cint
end

function WF_get_E_solution_length(self)
    @ccall LibLSE.WF_get_E_solution_length(self::Ptr{WaveFunction})::Cint
end

function WF_get_E_solution_data(self)
    @ccall LibLSE.WF_get_E_solution_data(self::Ptr{WaveFunction})::Ptr{Cdouble}
end

function psi_n(self, r, n, theta)
    @ccall LibLSE.psi_n(self::Ptr{WaveFunction}, r::Cdouble, n::Cint, theta::Cdouble)::Cint
end

function psi_n_ft(self, p, n)
    @ccall LibLSE.psi_n_ft(self::Ptr{WaveFunction}, p::Cdouble, n::Cint)::Cint
end

function psi_n_ftcomplex(self, complex)
    @ccall LibLSE.psi_n_ftcomplex(self::Ptr{WaveFunction}, complex::Cdouble)::Cint
end

function psi_n_batch(self, r_values, complex)
    @ccall LibLSE.psi_n_batch(self::Ptr{WaveFunction}, r_values::Ptr{Cdouble}, complex::Cdouble)::Cvoid
end

function psi_n_ft_batch(self, p_values, complex)
    @ccall LibLSE.psi_n_ft_batch(self::Ptr{WaveFunction}, p_values::Ptr{Cdouble}, complex::Cdouble)::Cvoid
end

function doublefactorial(n)
    @ccall LibLSE.doublefactorial(n::Cint)::Cdouble
end

function r_n(n)
    @ccall LibLSE.r_n(n::Cint)::Cdouble
end

function nu_n(n)
    @ccall LibLSE.nu_n(n::Cint)::Cdouble
end

function N_nl(n, l)
    @ccall LibLSE.N_nl(n::Cint, l::Cint)::Cdouble
end

function N_n_n_prime(n, n_prime, l)
    @ccall LibLSE.N_n_n_prime(n::Cint, n_prime::Cint, l::Cint)::Cdouble
end

function T_n_n_prime(n, n_prime, l)
    @ccall LibLSE.T_n_n_prime(n::Cint, n_prime::Cint, l::Cint)::Cdouble
end

function V_r_n_n_prime(n, n_prime, l)
    @ccall LibLSE.V_r_n_n_prime(n::Cint, n_prime::Cint, l::Cint)::Cdouble
end

function V_1_r_n_n_prime(n, n_prime, l)
    @ccall LibLSE.V_1_r_n_n_prime(n::Cint, n_prime::Cint, l::Cint)::Cdouble
end

function N_nl_tilde(n, l)
    @ccall LibLSE.N_nl_tilde(n::Cint, l::Cint)::Cdouble
end

function N_n_n_prime_tilde(n, n_prime, l)
    @ccall LibLSE.N_n_n_prime_tilde(n::Cint, n_prime::Cint, l::Cint)::Cdouble
end

function V_r_n_n_prime_tilde(n, n_prime, l)
    @ccall LibLSE.V_r_n_n_prime_tilde(n::Cint, n_prime::Cint, l::Cint)::Cdouble
end

function integrand_complex(r, complex)
    @ccall LibLSE.integrand_complex(r::Cdouble, complex::Cdouble)::Cint
end

function integrand(r, p, n, l)
    @ccall LibLSE.integrand(r::Cdouble, p::Cdouble, n::Cint, l::Cint)::Cint
end

const matrix = Cint

struct LSE
    Ngauss::Cint
    Lambda::Cdouble
    epsilon::Cdouble
    E::Cdouble
    T::Ptr{matrix}
    V::Ptr{matrix}
    G::Ptr{matrix}
    iIVG::Ptr{matrix}
    xi::Ptr{Cdouble}
    wi::Ptr{Cdouble}
    complex::Cdouble
    table::Ptr{Cint}
    wf::Ptr{WaveFunction}
    E_vec::Ptr{Cdouble}
end

function lse_malloc(pNgauss, Lambda, epsilon)
    @ccall LibLSE.lse_malloc(pNgauss::Cint, Lambda::Cdouble, epsilon::Cdouble)::Ptr{LSE}
end

function lse_compute(app, E)
    @ccall LibLSE.lse_compute(app::Ptr{LSE}, E::Cdouble)::Cint
end

function lse_free(app)
    @ccall LibLSE.lse_free(app::Ptr{LSE})::Cvoid
end

function lse_get_g_size(app, rows, cols)
    @ccall LibLSE.lse_get_g_size(app::Ptr{LSE}, rows::Ptr{Cuint}, cols::Ptr{Cuint})::Cvoid
end

function lse_get_v_size(app, rows, cols)
    @ccall LibLSE.lse_get_v_size(app::Ptr{LSE}, rows::Ptr{Cuint}, cols::Ptr{Cuint})::Cvoid
end

function lse_get_t_size(app, rows, cols)
    @ccall LibLSE.lse_get_t_size(app::Ptr{LSE}, rows::Ptr{Cuint}, cols::Ptr{Cuint})::Cvoid
end

function lse_get_ivg_size(app, rows, cols)
    @ccall LibLSE.lse_get_ivg_size(app::Ptr{LSE}, rows::Ptr{Cuint}, cols::Ptr{Cuint})::Cvoid
end

function lse_get_psi_size(self, rows, cols)
    @ccall LibLSE.lse_get_psi_size(self::Ptr{LSE}, rows::Ptr{Cuint}, cols::Ptr{Cuint})::Cvoid
end

function lse_gmat(self)
    @ccall LibLSE.lse_gmat(self::Ptr{LSE})::Cint
end

function lse_vmat(self)
    @ccall LibLSE.lse_vmat(self::Ptr{LSE})::Cint
end

function lse_tmat(self)
    @ccall LibLSE.lse_tmat(self::Ptr{LSE})::Cint
end

function lse_refresh(self, E)
    @ccall LibLSE.lse_refresh(self::Ptr{LSE}, E::Cdouble)::Cvoid
end

function xsqrt(x)
    @ccall LibLSE.xsqrt(x::Cint)::Cint
end

function square(x)
    @ccall LibLSE.square(x::Cdouble)::Cdouble
end

function V00(E, complex)
    @ccall LibLSE.V00(E::Cdouble, complex::Cdouble)::Cint
end

function V01(E, complex)
    @ccall LibLSE.V01(E::Cdouble, complex::Cdouble)::Cint
end

function V10(E, complex)
    @ccall LibLSE.V10(E::Cdouble, complex::Cdouble)::Cint
end

function V11(E, complex)
    @ccall LibLSE.V11(E::Cdouble, complex::Cdouble)::Cint
end

function V(alpha, beta, E, complex)
    @ccall LibLSE.V(alpha::Cint, beta::Cint, E::Cdouble, complex::Cdouble)::Cint
end

const PNGAUSS = 40

const RNGAUSS = 40

const RLAMBDA = 20

const f_pi = 0.092

const g_pi = 0.5704

const g_c = 1

const m_pi = 0.138039407

const m_K = 0.498

const m_eta = 0.548

const m_B = 5.27934

const m_B_star = 5.32471

const m_B_s = 5.36692

const m_B_star_s = 5.4154

const m11 = m_B

const m12 = m_B_star

const m21 = m_B_s

const m22 = m_B_star_s

const m_Xb13P = 10.5134 - (m11 + m12)

const m_Xb12P = 10.25546 - (m11 + m12)

const m_Xb11P = 9.89278 - (m11 + m12)

const a_l = 0.06426 * 5.068

const g = 0.898794378386677 ÷ a_l

const g1 = g * 0.014926616931653945

const g2 = g * 0.006467550544943349

const delta0 = 0

const delta1 = ((m21 + m22) - m11) - m12

const mu0 = (m11 * m12) ÷ (m11 + m12)

const mu1 = (m21 * m22) ÷ (m21 + m22)

const m_c = 1.85

const a_Cornell = 1.95

const g_qm = 1

const N_MAX = 30

const R_1 = 0.02 * 5.068

const R_N_MAX = 30 * 5.068

const C_T = ((-1 * 0.19732 * 0.19732) ÷ (0.5 * 4.18)) * 5.068

const PI = 3.141592653589793

const V0FIT = 0.4

const SIGMA_L = 0.0199179973550142

const SIGMA = (SIGMA_L ÷ a_l) ÷ a_l

const ALPHA = 0.476814032326273

# exports
const PREFIXES = ["lse_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
