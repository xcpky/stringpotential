import marimo

__generated_with = "0.12.8"
app = marimo.App(auto_download=["html"])


@app.cell
def _():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import figure
    import scipy
    from scipy import special
    from scipy.fft import fft,ifft,fftfreq
    from numba import njit
    return fft, fftfreq, figure, ifft, njit, np, plt, scipy, special


@app.cell
def _(plt):
    plt.rcParams["figure.figsize"] = (15,10)
    return


@app.cell
def _(math, np, special):
    n_max = 30
    #r_1 = 0.02 #fm
    #r_n_max = 30 #fm

    r_1 = 0.02*5.068 #GeV^-1
    r_n_max = 30*5.068 #GeV^-1

    #c_T = 2*41.47 #MeV fm^2 where the factor two comes from \mu=m/2 (in NN system)
    c_T = (-1)*0.19732**2/(1/2*4.18)*5.068
    #c_T = -2.0

    #define function for double factorial
    def doublefactorial(n):
        if n in (0, 1):
            return 1
        else:
            return n * doublefactorial(n-2)

    def r_n(n):
        return r_1*np.power(r_n_max/r_1,(n-1)/(n_max-1))

    def nu_n(n):
        return 1/(r_n(n)**2)

    def N_nl(n,l):
        return np.sqrt(2**(l+2)*(2*nu_n(n))**(l+3/2)/(np.sqrt(np.pi)*doublefactorial(2*l+1)))

    def N_n_n_prime(n,n_prime,l):
        return np.power(2*np.sqrt(nu_n(n)*nu_n(n_prime))/(nu_n(n)+nu_n(n_prime)),l+3/2)

    def T_n_n_prime(n,n_prime,l):
        return c_T*(2*l+3)*nu_n(n)*nu_n(n_prime)/(nu_n(n)+nu_n(n_prime))*N_n_n_prime(n,n_prime,l)

    def T_n_n_prime_test(n,n_prime,l):
        return (2*l*(l-1)/(2*l+1)*(nu_n(n)+nu_n(n_prime))-2*nu_n(n_prime)*(l+2)+2*nu_n(n_prime)**2*(2*l+3)/(nu_n(n)+nu_n(n_prime)))*N_n_n_prime(n,n_prime,l)

    def V_r_n_n_prime(n,n_prime,l):
        return 1/(np.sqrt(2*(nu_n(n)+nu_n(n_prime))))*doublefactorial(2*l+2)/doublefactorial(2*l+1)*N_n_n_prime(n,n_prime,l)

    def V_1_r_n_n_prime(n,n_prime,l):
        return 2/np.sqrt(np.pi)*np.power(2,l)*math.factorial(l)/(doublefactorial(2*l+1))*np.sqrt(nu_n(n)+nu_n(n_prime))*N_n_n_prime(n,n_prime,l)

    #for differnet normalization, i.e. psi->c*r^(l+1) for r->0
    #----------------------------
    def N_nl_tilde(n,l):
        return np.sqrt(2**(l+3)*(2*nu_n(n))**(l+5/2)/(np.sqrt(np.pi)*doublefactorial(2*l+3)))

    def N_n_n_prime_tilde(n,n_prime,l):
        return np.power(2*np.sqrt(nu_n(n)*nu_n(n_prime))/(nu_n(n)+nu_n(n_prime)),l+5/2)

    def T_n_n_prime_tilde_test(n,n_prime,l):
        return (2*((2*l+3)*nu_n(n)*nu_n(n_prime)-2*nu_n(n_prime)**2)/(nu_n(n)+nu_n(n_prime))-np.sqrt(2*(nu_n(n)+nu_n(n_prime)))*l*(l+1)*doublefactorial(2*l+2)/doublefactorial(2*l+3))*N_n_n_prime_tilde(n,n_prime,l)

    def T_n_n_prime_tilde(n,n_prime,l):
        return c_T*((2*l+5)*nu_n(n)*nu_n(n_prime)/(nu_n(n)+nu_n(n_prime))-(l+1))*N_n_n_prime_tilde(n,n_prime,l)

    def V_r_n_n_prime_tilde(n,n_prime,l):
        return 1/(np.sqrt(2*(nu_n(n)+nu_n(n_prime))))*doublefactorial(2*l+4)/doublefactorial(2*l+3)*N_n_n_prime_tilde(n,n_prime,l)

    def V_1_r_n_n_prime_tilde(n,n_prime,l):
        return doublefactorial(2*l+2)/(doublefactorial(2*l+3))*np.sqrt(2*(nu_n(n)+nu_n(n_prime)))*N_n_n_prime_tilde(n,n_prime,l)
    #----------------------------

    def V_e_mu_n_n_prime(n,n_prime,l,mu):
        return np.power(2*np.sqrt(nu_n(n)*nu_n(n_prime))/(nu_n(n)+nu_n(n_prime)+mu),l+3/2)

    def V_e_mu_lin_n_n_prime_l_0(n,n_prime,mu):
        return np.sqrt(4*(2*nu_n(n))**(3/2)/(np.sqrt(np.pi)))*np.sqrt(4*(2*nu_n(n_prime))**(3/2)/(np.sqrt(np.pi)))*(1/(2*(nu_n(n)+nu_n(n_prime)))-np.sqrt(np.pi)*mu*special.erfcx(mu/(2*np.sqrt(nu_n(n)+nu_n(n_prime))))/(4*(nu_n(n)+nu_n(n_prime))**(3/2)))

    def V_deuteron_n_n_prime(n,n_prime):
        return -626.885*V_e_mu_lin_n_n_prime_l_0(n,n_prime,1.55)+1438.72*V_e_mu_lin_n_n_prime_l_0(n,n_prime,3.11)
    return (
        N_n_n_prime,
        N_n_n_prime_tilde,
        N_nl,
        N_nl_tilde,
        T_n_n_prime,
        T_n_n_prime_test,
        T_n_n_prime_tilde,
        T_n_n_prime_tilde_test,
        V_1_r_n_n_prime,
        V_1_r_n_n_prime_tilde,
        V_deuteron_n_n_prime,
        V_e_mu_lin_n_n_prime_l_0,
        V_e_mu_n_n_prime,
        V_r_n_n_prime,
        V_r_n_n_prime_tilde,
        c_T,
        doublefactorial,
        n_max,
        nu_n,
        r_1,
        r_n,
        r_n_max,
    )


@app.cell
def _():
    #def constants
    f_pi = 0.092 #GeV
    g_pi = 0.5704 #GeV

    m_pi = 0.138039407 #GeV
    m_K = 0.498 #GeV
    m_eta = 0.548 #GeV

    m_B = 5.27934 #GeV
    m_B_star = 5.32471 #GeV
    m_B_s = 5.36692 #GeV
    m_B_star_s = 5.4154 #GeV

    m11 = m_B
    m12= m_B_star

    m21 = m_B_s
    m22 = m_B_star_s

    #constants to express Schroedinger equation in dimensionless units (according to Cornell paper)
    m_c=1.85 #GeV
    a_Cornell=1.95 #GeV^-1
    return (
        a_Cornell,
        f_pi,
        g_pi,
        m11,
        m12,
        m21,
        m22,
        m_B,
        m_B_s,
        m_B_star,
        m_B_star_s,
        m_K,
        m_c,
        m_eta,
        m_pi,
    )


@app.cell
def _(N_n_n_prime, T_n_n_prime, V_1_r_n_n_prime, V_r_n_n_prime, n_max, scipy):
    #define Schrödinger eq. as matrix eq.

    l = 0
    Lam = 0.0

    #scale factor in order to fit to \chi(b1)(1P)
    scale_factor=1.0

    #use data of D200 fit
    V_0_l = -0.08858455651140933 
    sigma_l = 0.0199179973550142
    alpha_l = 0.476814032326273
    r_0_l = 15

    a_l = 0.06426*5.068 #GeV^(-1)
    V_0 = V_0_l/a_l
    sigma = sigma_l/(a_l**2)
    alpha = alpha_l
    r_0 = r_0_l*a_l

    V_0_fit=0.4


    #H_deuteron = []
    H_Cornell = []
    N_n_n = []

    for i in range(n_max):
        n = i+1
    
        H_i = []
        N_n = []

        for j in range(n_max):
            n_prime = j+1

            #H_ij = T_n_n_prime(n,n_prime,0)+V_deuteron_n_n_prime(n,n_prime)
            #H_ij = T_n_n_prime(n,n_prime,l)+scale_factor*(V_0-sigma*r_0+alpha/r_0)*N_n_n_prime(n,n_prime,l)+sigma*V_r_n_n_prime(n,n_prime,l)-alpha*V_1_r_n_n_prime(n,n_prime,l)
            H_ij = -T_n_n_prime(n,n_prime,l)+V_0_fit*N_n_n_prime(n,n_prime,l)+sigma*V_r_n_n_prime(n,n_prime,l)-alpha*V_1_r_n_n_prime(n,n_prime,l)
            #H_ij = -T_n_n_prime(n,n_prime,l)+V_r_n_n_prime(n,n_prime,l)-Lam*V_1_r_n_n_prime(n,n_prime,l)
            N_nn = N_n_n_prime(n,n_prime,l)
        
            H_i.append(H_ij)
            N_n.append(N_nn)
    
        #H_deuteron.append(H_i)
        H_Cornell.append(H_i)
        N_n_n.append(N_n)

    #use scipy to solve generalised eigenvalue problem
    E_solution, c_solution = scipy.linalg.eigh(H_Cornell,b=N_n_n)
    return (
        E_solution,
        H_Cornell,
        H_i,
        H_ij,
        Lam,
        N_n,
        N_n_n,
        N_nn,
        V_0,
        V_0_fit,
        V_0_l,
        a_l,
        alpha,
        alpha_l,
        c_solution,
        i,
        j,
        l,
        n,
        n_prime,
        r_0,
        r_0_l,
        scale_factor,
        sigma,
        sigma_l,
    )


@app.cell
def _(V_0, alpha, r_0, sigma):
    print(V_0-sigma*r_0+alpha/r_0)
    return


@app.cell
def _(E_solution, np):
    np.savetxt('E_solution_GaussBasis_Cornell_30_V0_04.txt',E_solution,fmt='%.8e')
    return


@app.cell
def _(c_solution, np):
    np.savetxt('c_solution_GaussBasis_Cornell_30_V0_04.txt',c_solution,fmt='%.8e')
    return


@app.cell
def _(m11, m12):
    m_Xb11P=9.89278-(m11+m12)
    print(m_Xb11P)
    return (m_Xb11P,)


@app.cell
def _(E_solution):
    print(E_solution)
    return


@app.cell
def _(E_solution, a_Cornell, m_c):
    print('In dimensionless units, \zeta:',m_c**(-1)*(m_c*a_Cornell)**(4/3)*E_solution)
    return


@app.cell
def _(c_solution):
    print(c_solution[:,0])
    return


@app.cell
def _(N_nl, N_nl_tilde, c_solution, n_max, np, nu_n, special):
    def psi_nlm(r, n, l, m, theta, phi):
        psi = 0.0
        for _i in range(n_max):
            _n_prime = _i + 1
            psi = psi + c_solution[:, n][_i] * N_nl(_n_prime, _l) * r ** _l * np.exp(-nu_n(_n_prime) * r ** 2) * special.sph_harm(m, _l, theta, phi)
        return psi

    def psi_nlm_tilde(r, n, l, m, theta, phi):
        psi = 0
        for _i in range(n_max):
            _n_prime = _i + 1
            psi = psi + c_solution[:, n][_i] * N_nl_tilde(_n_prime, _l) * r ** (_l + 1) * np.exp(-nu_n(_n_prime) * r ** 2) * special.sph_harm(m, _l, theta, phi)
        return psi

    def psi_nlm_cs(r, n, l, m, theta, phi):
        if psi_nlm(0.0 + 0.0001, n, _l, m, theta, phi) < 0.0:
            return -psi_nlm(r, n, _l, m, theta, phi)
        else:
            return psi_nlm(r, n, _l, m, theta, phi)
    return psi_nlm, psi_nlm_cs, psi_nlm_tilde


@app.cell
def _(E_solution, c_solution):
    res = 0.0
    for _i in range(E_solution.size):
        res = res + c_solution[:, 0][_i] ** 2
    return (res,)


@app.cell
def _(res):
    print(res)
    return


@app.cell
def _(c_solution, n_max, np, nu_n, psi_nlm):
    def N_nl_p(n_prime):
        return np.sqrt(4 * np.sqrt(2) / (9 * np.pi ** (7 / 2)) * nu_n(_n_prime) ** (5 / 2))

    def psi_nlm_p(p, n):
        psi_p = 0.0 + 0j
        for _i in range(n_max):
            _n_prime = _i + 1
            psi_p = psi_p + c_solution[:, n][_i] * N_nl_p(_n_prime) * np.sqrt(3) * np.pi * 1j * 1 / (4 * nu_n(_n_prime) ** (5 / 2)) * p * np.exp(-p ** 2 / (4 * nu_n(_n_prime)))
        if psi_nlm(0.0 + 0.0001, n, 1, 0, 0, 0) < 0.0:
            return -psi_p
        else:
            return psi_p
    return N_nl_p, psi_nlm_p


@app.cell
def _(np, plt, psi_nlm_cs, special):
    r_plot = np.arange(0.0, 7.0, 0.0001)
    _l = 0
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 0, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 1, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 2, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 3, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 4, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 5, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot, r_plot * psi_nlm_cs(r_plot, 6, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.grid()
    return (r_plot,)


@app.cell
def _(np):
    Cornell_data_1S_x=np.loadtxt('Cornell_data_1S.txt',usecols=0)
    Cornell_data_1S_y=np.loadtxt('Cornell_data_1S.txt',usecols=1)

    Cornell_data_2S_x=np.loadtxt('Cornell_data_2S.txt',usecols=0)
    Cornell_data_2S_y=np.loadtxt('Cornell_data_2S.txt',usecols=1)

    Cornell_data_3S_x=np.loadtxt('Cornell_data_3S.txt',usecols=0)
    Cornell_data_3S_y=np.loadtxt('Cornell_data_3S.txt',usecols=1)

    Cornell_data_4S_x=np.loadtxt('Cornell_data_4S.txt',usecols=0)
    Cornell_data_4S_y=np.loadtxt('Cornell_data_4S.txt',usecols=1)
    return (
        Cornell_data_1S_x,
        Cornell_data_1S_y,
        Cornell_data_2S_x,
        Cornell_data_2S_y,
        Cornell_data_3S_x,
        Cornell_data_3S_y,
        Cornell_data_4S_x,
        Cornell_data_4S_y,
    )


@app.cell
def _():
    import scipy.interpolate as sciint
    return (sciint,)


@app.cell
def _(
    Cornell_data_1S_x,
    Cornell_data_1S_y,
    Cornell_data_2S_x,
    Cornell_data_2S_y,
    Cornell_data_3S_x,
    Cornell_data_3S_y,
    Cornell_data_4S_x,
    Cornell_data_4S_y,
    sciint,
):
    test_func_1S=sciint.CubicSpline(Cornell_data_1S_x, Cornell_data_1S_y)
    test_func_2S=sciint.CubicSpline(Cornell_data_2S_x, Cornell_data_2S_y)
    test_func_3S=sciint.CubicSpline(Cornell_data_3S_x, Cornell_data_3S_y)
    test_func_4S=sciint.CubicSpline(Cornell_data_4S_x, Cornell_data_4S_y)
    return test_func_1S, test_func_2S, test_func_3S, test_func_4S


@app.cell
def _(
    a_Cornell,
    m_c,
    np,
    plt,
    psi_nlm_cs,
    special,
    test_func_1S,
    test_func_2S,
    test_func_3S,
    test_func_4S,
):
    x_interpolate=np.linspace(0.0,6.0,1000)
    rho=a_Cornell**(-1)*(m_c*a_Cornell)**(1/3)*x_interpolate

    yu0=x_interpolate*psi_nlm_cs(x_interpolate,0,0,0,0,0).real/special.sph_harm(0,0,0,0).real+0.02
    yd0=x_interpolate*psi_nlm_cs(x_interpolate,0,0,0,0,0).real/special.sph_harm(0,0,0,0).real-0.02

    yu1=x_interpolate*psi_nlm_cs(x_interpolate,1,0,0,0,0).real/special.sph_harm(0,0,0,0).real+0.02
    yd1=x_interpolate*psi_nlm_cs(x_interpolate,1,0,0,0,0).real/special.sph_harm(0,0,0,0).real-0.02

    yu2=x_interpolate*psi_nlm_cs(x_interpolate,2,0,0,0,0).real/special.sph_harm(0,0,0,0).real+0.02
    yd2=x_interpolate*psi_nlm_cs(x_interpolate,2,0,0,0,0).real/special.sph_harm(0,0,0,0).real-0.02

    yu3=x_interpolate*psi_nlm_cs(x_interpolate,3,0,0,0,0).real/special.sph_harm(0,0,0,0).real+0.02
    yd3=x_interpolate*psi_nlm_cs(x_interpolate,3,0,0,0,0).real/special.sph_harm(0,0,0,0).real-0.02

    x_interpolate_scale=0.94

    plt.plot(x_interpolate,x_interpolate*psi_nlm_cs(x_interpolate,0,0,0,0,0).real/special.sph_harm(0,0,0,0).real,color='royalblue',linestyle='dashed')
    plt.plot(x_interpolate,x_interpolate*psi_nlm_cs(x_interpolate,1,0,0,0,0).real/special.sph_harm(0,0,0,0).real,color='orange',linestyle='dashed')
    plt.plot(x_interpolate,x_interpolate*psi_nlm_cs(x_interpolate,2,0,0,0,0).real/special.sph_harm(0,0,0,0).real,color='green',linestyle='dashed')
    plt.plot(x_interpolate,x_interpolate*psi_nlm_cs(x_interpolate,3,0,0,0,0).real/special.sph_harm(0,0,0,0).real,color='firebrick',linestyle='dashed')

    #plt.fill_between(x_interpolate,y1=yu0,y2=yd0,color='blue',alpha=0.2)
    #plt.fill_between(x_interpolate,y1=yu1,y2=yd1,color='orange',alpha=0.2)
    #plt.fill_between(x_interpolate,y1=yu2,y2=yd2,color='green',alpha=0.2)
    #plt.fill_between(x_interpolate,y1=yu3,y2=yd3,color='firebrick',alpha=0.2)

    plt.plot(x_interpolate,test_func_1S(x_interpolate),color='royalblue')
    plt.plot(x_interpolate,test_func_2S(x_interpolate),color='orange')
    plt.plot(x_interpolate,test_func_3S(x_interpolate),color='green')
    plt.plot(x_interpolate,test_func_4S(x_interpolate),color='firebrick')

    plt.text(2.0,0.68,r'$1S$',fontsize=20,color='royalblue')
    plt.text(1.3,0.4,r'$2S$',fontsize=20,color='orange')
    plt.text(4.5,0.65,r'$3S$',fontsize=20,color='green')
    plt.text(4.3,0.2,r'$4S$',fontsize=20,color='firebrick')

    plt.xlabel(r'$\rho$',fontsize=20)
    plt.ylabel(r'$u(\rho)$',fontsize=20)
    plt.grid()
    plt.savefig('Cornell_comp.pdf',dpi=500)
    return (
        rho,
        x_interpolate,
        x_interpolate_scale,
        yd0,
        yd1,
        yd2,
        yd3,
        yu0,
        yu1,
        yu2,
        yu3,
    )


@app.cell
def _(
    Cornell_data_1S_x,
    Cornell_data_1S_y,
    Cornell_data_2S_x,
    Cornell_data_2S_y,
    Cornell_data_3S_x,
    Cornell_data_3S_y,
    Cornell_data_4S_x,
    Cornell_data_4S_y,
    plt,
):
    plt.plot(Cornell_data_1S_x,Cornell_data_1S_y)
    plt.plot(Cornell_data_2S_x,Cornell_data_2S_y)
    plt.plot(Cornell_data_3S_x,Cornell_data_3S_y)
    plt.plot(Cornell_data_4S_x,Cornell_data_4S_y)
    return


@app.cell
def _(np, plt, psi_nlm_cs, special):
    r_plot_1 = np.arange(0.0, 7.0, 0.0001)
    _l = 0
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 0, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 1, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 2, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 3, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 4, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 5, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.plot(r_plot_1, r_plot_1 * psi_nlm_cs(r_plot_1, 6, _l, 0, 0, 0).real / special.sph_harm(0, 0, 0, 0).real)
    plt.grid()
    return (r_plot_1,)


@app.cell
def _(np):
    def PIgauss(x, a1, a2, mu1, mu2, sigma1, sigma2):
        return a1 * np.exp(-(_x - mu1) ** 2 / (2 * sigma1 ** 2)) + a2 * np.exp(-(_x - mu2) ** 2 / (2 * sigma2 ** 2))
    return (PIgauss,)


@app.cell
def _(PIgauss, plt, r_plot_1):
    plt.plot(r_plot_1, PIgauss(r_plot_1, 1, -1, 1, 4, 0.5, 2))
    return


@app.cell
def _(E_solution, plt, psi_nlm, r_plot_1):
    _fig, _axs = plt.subplots(1, 4, figsize=(15, 15), sharey=True)
    _axs[0].plot(r_plot_1, psi_nlm(r_plot_1, 0, 0, 0, 0, 0))
    _axs[0].text(4, 0.19, '$E_0=%.2f \\,\\mathrm{GeV}$' % E_solution[0], fontsize=15)
    _axs[0].set_ylabel('$\\phi_{i,l=0}(r)\\,/\\,\\mathrm{GeV}^{3/2}$', fontsize=20)
    _axs[0].set_xlabel('$r\\,/\\,\\mathrm{GeV}^{-1}$', fontsize=20)
    _axs[0].grid()
    _axs[1].plot(r_plot_1, psi_nlm(r_plot_1, 1, 0, 0, 0, 0))
    _axs[1].text(4, 0.19, '$E_1=%.2f \\,\\mathrm{GeV}$' % E_solution[1], fontsize=15)
    _axs[1].set_xlabel('$r\\,/\\,\\mathrm{GeV}^{-1}$', fontsize=20)
    _axs[1].grid()
    _axs[2].plot(r_plot_1, psi_nlm(r_plot_1, 2, 0, 0, 0, 0))
    _axs[2].text(4, 0.19, '$E_2=%.2f \\,\\mathrm{GeV}$' % E_solution[2], fontsize=15)
    _axs[2].set_xlabel('$r\\,/\\,\\mathrm{GeV}^{-1}$', fontsize=20)
    _axs[2].grid()
    _axs[3].plot(r_plot_1, psi_nlm(r_plot_1, 3, 0, 0, 0, 0))
    _axs[3].text(5, 0.19, '$E_3=%.2f \\,\\mathrm{GeV}$' % E_solution[3], fontsize=15)
    _axs[3].set_xlabel('$r\\,/\\,\\mathrm{GeV}^{-1}$', fontsize=20)
    _axs[3].grid()
    _fig.subplots_adjust(wspace=0, hspace=0)
    return


@app.cell
def _(np):
    def Gauss_mesh(n):
        """
        return xi(n),wi(n)
        """
        out = np.polynomial.legendre.leggauss(n)
        return (out[0], out[1])

    def Gauss_Leg_mapping(x, a, b):
        return (b - a) * _x / 2 + (a + b) / 2

    def Gauss_Legendre(f, a, b, x, w, *args):
        """
        f=function(x_int,*args)
        arg=f,a,b,x,w,*args
        """
        res = 0.0 + 0j
        x_var = Gauss_Leg_mapping(_x, a, b)
        for _i in range(x_var.size):
            res = res + _w[_i] * f(x_var[_i], *args)
        return (b - a) * res / 2
    return Gauss_Leg_mapping, Gauss_Legendre, Gauss_mesh


@app.cell
def _(Gauss_Legendre, Gauss_mesh, N_nl_p, np, nu_n):
    _x, _w = Gauss_mesh(80)
    _n_prime = 3

    def integrandTestpsip(p):
        return 4 * np.pi * p ** 2 * np.abs(N_nl_p(_n_prime) * np.sqrt(3) * np.pi * 1j * 1 / (4 * nu_n(_n_prime) ** (5 / 2)) * p * np.exp(-p ** 2 / (4 * nu_n(_n_prime)))) ** 2
    print(Gauss_Legendre(integrandTestpsip, 0, 1000.0, _x, _w))
    return (integrandTestpsip,)


@app.cell
def _(m11, m12, m21, m22, np):
    mu = np.array([m11*m12/(m11+m12),m21*m22/(m21+m22)])

    Delta = np.array([m11+m12,m21+m22])-(m11+m12) #m_1^0=m11, m_2^0=m12

    def q_0_func(E):
        q0 = np.array([np.sqrt(2*mu[0]*(E-Delta[0])+0j),np.sqrt(2*mu[1]*(E-Delta[1])+0j)])
        return q0.real+0j
    return Delta, mu, q_0_func


@app.cell
def _(np):
    def Gauss_herm_mesh(n):
        out = np.polynomial.hermite.hermgauss(n)
        return (out[0], out[1])

    def Gauss_herm(f, x, w, *args):
        res = 0.0 + 0j
        for _i in range(_x.size):
            res = res + _w[_i] * f(np.abs(_x[_i]), *args)
        return res / 2
    return Gauss_herm, Gauss_herm_mesh


@app.cell
def _(Gauss_Legendre, Gauss_herm, Gauss_herm_mesh, Gauss_mesh, np, psi_nlm_cs):
    def integrand_3S1(x, p, n):
        if p == 0:
            p = p + 1e-10
        else:
            p = p
        return np.exp(1j * _x * p) * _x / p * psi_nlm_cs(_x, n, 0, 0, 0, 0)

    def integrand_3P1(x, p, n):
        if p == 0:
            p = p + 1e-10
        else:
            p = p
        return np.exp(1j * _x * p) * _x / p * psi_nlm_cs(_x, n, 1, 0, 0, 0)

    def Ftpsi_3P1(p, n, N_gauss, Lambda):
        _x, _w = Gauss_mesh(N_gauss)
        result1 = -1j * (2 * np.pi) * (Gauss_Legendre(integrand_3P1, 0, Lambda, _x, _w, p, n) + Gauss_Legendre(integrand_3P1, 0, Lambda, _x, _w, -p, n))
        result2 = result1.real ** 2 + result1.imag ** 2
        return (result1, result2)

    def Ftpsi_S_1_3(p, n, N_gauss, Lambda):
        _x, _w = Gauss_mesh(N_gauss)
        result1 = -1j * (2 * np.pi) * (Gauss_Legendre(integrand_3S1, 0, Lambda, _x, _w, p, n) + Gauss_Legendre(integrand_3S1, 0, Lambda, _x, _w, -p, n))
        result2 = result1.real ** 2 + result1.imag ** 2
        return (result1, result2)

    def Ftpsi_S_1_3_hermite(p, n, N_gauss):
        _x, _w = Gauss_herm_mesh(N_gauss)
        result1 = -1j * (2 * np.pi) * (Gauss_herm(integrand_3S1, _x, _w, p, n) + Gauss_herm(integrand_3S1, _x, _w, -p, n))
        result2 = result1.real ** 2 + result1.imag ** 2
        return (result1, result2)
    return (
        Ftpsi_3P1,
        Ftpsi_S_1_3,
        Ftpsi_S_1_3_hermite,
        integrand_3P1,
        integrand_3S1,
    )


@app.cell
def progressbar():
    def progressbar(i_in, imax, state):
        """
        state=np.array([False,False,...])
        """
        r = int(round(i_in * 100 / imax))
        for _i in range(state.size):
            if state[_i] == False and r >= round((_i + 1) * 100 / state.size):
                print(str(_i + 1) + '/' + str(state.size))
                state[_i] = True
        return state
    return (progressbar,)


@app.cell
def _(Delta, Gauss_mesh, np):
    E_array = np.linspace(Delta[0] - 1.0, Delta[1] + 0.4, 1000)
    print(len(E_array))
    Nq = 40
    Ntower = 30
    xq, _wq = Gauss_mesh(Nq)
    psi_mat_q1 = np.zeros((Ntower * len(E_array), Nq + 1))
    psi_mat_q2 = np.zeros((Ntower * len(E_array), Nq + 1))
    state = np.zeros(Ntower)
    return E_array, Nq, Ntower, psi_mat_q1, psi_mat_q2, state, xq


@app.cell
def _(
    E_array,
    Ftpsi_S_1_3,
    Nq,
    Ntower,
    progressbar,
    psi_mat_q1,
    psi_mat_q2,
    q_0_func,
    state,
    xq,
):
    for _i in range(Ntower):
        state_1 = progressbar(_i, Ntower, state)
        for _k in range(psi_mat_q1.shape[1] - 1):
            psi_mat_q1[_i, _k] = Ftpsi_S_1_3(xq[_k], _i, Nq, 20.0)[0].real
            psi_mat_q2[_i, _k] = psi_mat_q1[_i, _k]
        for _l in range(len(E_array)):
            psi_mat_q1[_i + Ntower * _l, Nq] = Ftpsi_S_1_3(q_0_func(E_array[_l])[0], _i, Nq, 20.0)[0].real
            psi_mat_q2[_i + Ntower * _l, Nq] = Ftpsi_S_1_3(q_0_func(E_array[_l])[1], _i, Nq, 20.0)[0].real
            psi_mat_q1[_i + Ntower * _l, :-1] = psi_mat_q1[_i, :-1]
            psi_mat_q2[_i + Ntower * _l, :-1] = psi_mat_q2[_i, :-1]
    return (state_1,)


@app.cell
def _(
    E_array,
    Nq,
    Ntower,
    progressbar,
    psi_mat_q1,
    psi_mat_q2,
    psi_nlm_p,
    q_0_func,
    state_1,
    xq,
):
    for _i in range(Ntower):
        state_2 = progressbar(_i, Ntower, state_1)
        for _k in range(psi_mat_q1.shape[1] - 1):
            psi_mat_q1[_i, _k] = psi_nlm_p(xq[_k], _i).real
            psi_mat_q2[_i, _k] = psi_mat_q1[_i, _k]
        for _l in range(len(E_array)):
            psi_mat_q1[_i + Ntower * _l, Nq] = psi_nlm_p(q_0_func(E_array[_l])[0], _i).real
            psi_mat_q2[_i + Ntower * _l, Nq] = psi_nlm_p(q_0_func(E_array[_l])[1], _i).real
            psi_mat_q1[_i + Ntower * _l, :-1] = psi_mat_q1[_i, :-1]
            psi_mat_q2[_i + Ntower * _l, :-1] = psi_mat_q2[_i, :-1]
    return (state_2,)


@app.cell
def _(np, psi_mat_q1, psi_mat_q2):
    np.savetxt('psi_mat_q1_Nq40_E1000_104_Ntower30.txt',psi_mat_q1,fmt='%.8e')
    np.savetxt('psi_mat_q2_Nq40_E1000_104_Ntower30.txt',psi_mat_q2,fmt='%.8e')
    return


@app.cell
def _(psi_mat_q1):
    print(psi_mat_q1)
    #print(psi_mat_q2[29000,40])
    return


@app.cell
def _(psi_mat_q1):
    psi_mat_q1.shape
    return


@app.cell
def _(plt, psi_mat_q1, xq):
    plt.plot(xq,psi_mat_q1[0,:-1])
    plt.plot(xq,psi_mat_q1[1,:-1])
    plt.plot(xq,psi_mat_q1[2,:-1])
    plt.plot(xq,psi_mat_q1[3,:-1])

    plt.plot(xq,psi_mat_q1[4,:-1])
    plt.plot(xq,psi_mat_q1[5,:-1])
    plt.plot(xq,psi_mat_q1[6,:-1])
    plt.plot(xq,psi_mat_q1[7,:-1])

    plt.plot(xq,psi_mat_q1[8,:-1])
    plt.plot(xq,psi_mat_q1[9,:-1])
    plt.plot(xq,psi_mat_q1[10,:-1])
    plt.plot(xq,psi_mat_q1[11,:-1])
    #plt.xlim(-0.1,0.0)
    return


@app.cell
def _(Ftpsi_S_1_3, np, q_0_func):
    E_array_1 = np.linspace(-1.0, 1.0, 100)
    plot_psi_E_q1 = np.zeros(E_array_1.size)
    plot_psi_E_q2 = np.zeros(E_array_1.size)
    for _l in range(E_array_1.size):
        plot_psi_E_q1[_l] = Ftpsi_S_1_3(q_0_func(E_array_1[_l])[0], 0, 40, 20.0)
        plot_psi_E_q2[_l] = Ftpsi_S_1_3(q_0_func(E_array_1[_l])[1], 0, 40, 20.0)
    return E_array_1, plot_psi_E_q1, plot_psi_E_q2


@app.cell
def _(E_array_1, plot_psi_E_q1, plot_psi_E_q2, plt):
    plt.plot(E_array_1, plot_psi_E_q1)
    plt.plot(E_array_1, plot_psi_E_q2)
    return


@app.cell
def rectangular_int():
    def rectangular_int(f, a, b, N, *args):
        _w = (b - a) / N
        res = 0.0 + 0j
        for _i in range(N):
            x_i = a + _w * _i
            res = res + _w * f(x_i, *args)
        return res
    return (rectangular_int,)


@app.cell
def _(
    Gauss_Legendre,
    Gauss_herm,
    Gauss_herm_mesh,
    Gauss_mesh,
    rectangular_int,
):
    def test(x):
        return _x
    print(rectangular_int(test, 0, 10, 20))
    _x, _w = Gauss_mesh(20)
    print(Gauss_Legendre(test, 0, 10, _x, _w))
    _x, _w = Gauss_herm_mesh(100)
    print(Gauss_herm(test, _x, _w))
    return (test,)


@app.cell
def _(integrand, np, rectangular_int):
    def Ftpsi_S_1_3_rec(p,n,N,Lambda):
        result1 = -1j*(2*np.pi)*(rectangular_int(integrand,0,Lambda,N,p,n)+rectangular_int(integrand,0,Lambda,N,-p,n))
        result2 = result1.real**2+result1.imag**2
        return result1,result2
    return (Ftpsi_S_1_3_rec,)


@app.cell
def _(Ftpsi_S_1_3, Gauss_mesh, np):
    plot_p = np.linspace(-10.0, 10.0, 1000) + 0j
    plot_E = np.linspace(-2.5, 10.0, 1000)
    xq_1, _wq = Gauss_mesh(40)
    plot_psi_p_40_0 = np.zeros(plot_p.size, dtype=np.complex128)
    plot_psi_p_40_1 = np.zeros(plot_p.size, dtype=np.complex128)
    plot_psi_p_40_2 = np.zeros(plot_p.size, dtype=np.complex128)
    for _i in range(plot_p.size):
        plot_psi_p_40_0[_i] = Ftpsi_S_1_3(plot_p[_i], 0, 40, 20.0)[0]
        plot_psi_p_40_1[_i] = Ftpsi_S_1_3(plot_p[_i], 1, 40, 20.0)[0]
        plot_psi_p_40_2[_i] = Ftpsi_S_1_3(plot_p[_i], 2, 40, 20.0)[0]
    return (
        plot_E,
        plot_p,
        plot_psi_p_40_0,
        plot_psi_p_40_1,
        plot_psi_p_40_2,
        xq_1,
    )


@app.cell
def _(plot_psi_p_40_real_herm_xq, plt, xq_1):
    plt.plot(xq_1, plot_psi_p_40_real_herm_xq)
    return


@app.cell
def _(np):
    plot_p_1 = np.linspace(-10, 10, 1000)
    plot_p60 = np.linspace(-10, 10, 100)
    return plot_p60, plot_p_1


@app.cell
def _(N_nl, N_nl_p, np, plot_p_1, plot_psi_p_40_0):
    N_test = np.zeros(plot_p_1.size)
    for _i in range(plot_p_1.size):
        N_test[_i] = plot_p_1[_i] ** 2 * np.abs(N_nl_p(0) / N_nl(0, 0) * plot_psi_p_40_0[_i]) ** 2 * 20 / 1000
    print(np.sum(N_test))
    return (N_test,)


@app.cell
def _(N_nl, N_nl_p):
    print(N_nl_p(10)/N_nl(10,0)-N_nl_p(30)/N_nl(30,0))
    return


@app.cell
def _(np, plot_p_1, plt, psi_nlm_cs, psi_nlm_p):
    _fig, _axs = plt.subplots(1, 2, figsize=(15, 10))
    _axs[0].plot(plot_p_1, np.abs(psi_nlm_p(plot_p_1, 0)) ** 2, label='$|\\psi_{0}|^2$, analytical', color='darkblue', linestyle='dashed')
    _axs[0].plot(plot_p_1, np.abs(psi_nlm_p(plot_p_1, 1)) ** 2, label='$|\\psi_{1}|^2$, analytical', color='darkorange', linestyle='dashed')
    _axs[0].plot(plot_p_1, np.abs(psi_nlm_p(plot_p_1, 2)) ** 2, label='$|\\psi_{2}|^2$, analytical', color='darkgreen', linestyle='dashed')
    _axs[0].set_ylabel('$\\psi_{i,l=1}(p)\\,/\\,\\mathrm{GeV}^{-3/2}$', fontsize=20)
    _axs[0].set_xlabel('$p\\,/\\,\\mathrm{GeV}$', fontsize=20)
    _axs[0].set_xlim(-10.0, 10.0)
    _axs[0].legend()
    _axs[0].grid()
    plot_x = np.linspace(-2.0, 2.0, 1000)
    _axs[1].plot(plot_x, np.abs(psi_nlm_cs(plot_x, 0, 1, 0, 0, 0)) ** 2)
    _axs[1].plot(plot_x, np.abs(psi_nlm_cs(plot_x, 1, 1, 0, 0, 0)) ** 2)
    _axs[1].plot(plot_x, np.abs(psi_nlm_cs(plot_x, 2, 1, 0, 0, 0)) ** 2)
    _axs[1].set_ylabel('$\\psi_{i,l=1}(x)\\,/\\,\\mathrm{GeV}^{3/2}$', fontsize=20)
    _axs[1].set_xlabel('$x\\,/\\,\\mathrm{GeV}^{-1}$', fontsize=20)
    _axs[1].grid()
    return (plot_x,)


@app.cell
def _(Gauss_mesh, np, plt, psi_nlm_p):
    xq_2, _wq = Gauss_mesh(40)
    plot_p_Test = np.zeros(xq_2.size, dtype=np.complex128)
    for _i in range(xq_2.size):
        plot_p_Test[_i] = np.abs(psi_nlm_p(xq_2[_i], 0)) ** 2
    plt.plot(xq_2, plot_p_Test)
    return plot_p_Test, xq_2


@app.cell
def _(Ftpsi_S_1_3, np, plot_p_1):
    plot_psi_p_0_abs = np.zeros(len(plot_p_1), dtype=np.complex128)
    plot_psi_p_1_abs = np.zeros(len(plot_p_1), dtype=np.complex128)
    plot_psi_p_2_abs = np.zeros(len(plot_p_1), dtype=np.complex128)
    for _i in range(len(plot_p_1)):
        plot_psi_p_0_abs[_i] = Ftpsi_S_1_3(plot_p_1[_i], 0, 40, 20.0)[1]
        plot_psi_p_1_abs[_i] = Ftpsi_S_1_3(plot_p_1[_i], 1, 40, 20.0)[1]
        plot_psi_p_2_abs[_i] = Ftpsi_S_1_3(plot_p_1[_i], 2, 40, 20.0)[1]
    return plot_psi_p_0_abs, plot_psi_p_1_abs, plot_psi_p_2_abs


@app.cell
def _(Ftpsi_S_1_3, np, plot_p_1):
    plot_psi_p_3_abs = np.zeros(len(plot_p_1), dtype=np.complex128)
    plot_psi_p_4_abs = np.zeros(len(plot_p_1), dtype=np.complex128)
    plot_psi_p_5_abs = np.zeros(len(plot_p_1), dtype=np.complex128)
    for _i in range(len(plot_p_1)):
        plot_psi_p_3_abs[_i] = Ftpsi_S_1_3(plot_p_1[_i], 3, 40, 20.0)[1]
        plot_psi_p_4_abs[_i] = Ftpsi_S_1_3(plot_p_1[_i], 4, 40, 20.0)[1]
        plot_psi_p_5_abs[_i] = Ftpsi_S_1_3(plot_p_1[_i], 5, 40, 20.0)[1]
    return plot_psi_p_3_abs, plot_psi_p_4_abs, plot_psi_p_5_abs


@app.cell
def _(
    plot_p_1,
    plot_psi_p_0_abs,
    plot_psi_p_1_abs,
    plot_psi_p_2_abs,
    plot_psi_p_3_abs,
    plot_psi_p_4_abs,
    plot_psi_p_5_abs,
    plt,
):
    plt.plot(plot_p_1, plot_psi_p_0_abs.real / plot_psi_p_0_abs.real[499])
    plt.plot(plot_p_1, plot_psi_p_1_abs.real / plot_psi_p_1_abs.real[499])
    plt.plot(plot_p_1, plot_psi_p_2_abs.real / plot_psi_p_2_abs.real[499])
    plt.plot(plot_p_1, plot_psi_p_3_abs.real / plot_psi_p_3_abs.real[499])
    plt.plot(plot_p_1, plot_psi_p_4_abs.real / plot_psi_p_4_abs.real[499])
    plt.plot(plot_p_1, plot_psi_p_5_abs.real / plot_psi_p_5_abs.real[499])
    return


@app.cell
def _(np):
    m = 0.138039407 #GeV
    def test_q(q):
        return 1/(q**2+m**2)

    def test_r(r):
        return 1/(4*np.pi)*np.exp(-r*m)/r
    return m, test_q, test_r


@app.cell
def _(np, plt, test_r):
    _x = np.linspace(-10, 10, 1000)
    plt.plot(_x, test_r(_x))
    return


@app.cell
def _(Gauss_Legendre, Gauss_mesh, np, test_r):
    p = np.linspace(-1.0, 1.0, 1000)

    def integrand2(r, p):
        if p == 0:
            p = p + 1e-10
        else:
            p = p
        return np.exp(1j * r * p) * r / p * test_r(r)

    def integrand3(r, p):
        return np.exp(1j * r * p) * r / p * np.exp(-r ** 2)

    def ftTest(p, f, N_gauss, Lambda):
        _x, _w = Gauss_mesh(N_gauss)
        result1 = -1j * (Gauss_Legendre(f, 0, Lambda, _x, _w, p) + Gauss_Legendre(f, 0, Lambda, _x, _w, -p))
        result2 = abs(result1)
        return (result1, result2)
    return ftTest, integrand2, integrand3, p


@app.cell
def _(ftTest, integrand2, np, p):
    test_p = []
    for _i in range(len(p)):
        test_p.append(2 * np.pi * ftTest(p[_i], integrand2, 40, 50.0)[0].real)
    return (test_p,)


@app.cell
def _(p, plt, test_p, test_q):
    plt.plot(p,test_p)
    plt.plot(p,test_q(p),color='firebrick',linestyle='dashed')
    plt.grid()
    return


if __name__ == "__main__":
    app.run()
