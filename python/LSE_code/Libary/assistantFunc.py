import numpy as np
import matplotlib.pyplot as plt
from scipy import special

from Libary.constants import *

#define functions for Gauss-Legendre integration
#------------------------------------------------
def Gauss_mesh(n):
    '''
    return xi(n),wi(n)
    '''
    out=np.polynomial.legendre.leggauss(n)
    return out[0],out[1]

@njit
def Gauss_Leg_mapping(x,w,a,b):
    out1=(b-a)*x/2+(a+b)/2
    out2=(b-a)/2*w
    return out1,out2
    
@njit
def Gauss_Legendre(f,a,b,x,w,*args):
    '''
    f=function(x_int,*args)
    '''
    res=0.+0j
    x_var,w=Gauss_Leg_mapping(x,w,a,b)
    for i in range(x_var.size):
        res+=w[i]*f(x_var[i],*args)
    return res
#------------------------------------------------
#define Gauss points and weights
#------------------------------------------------

Nq_split=int(Nq/2)
qGauss,wGauss=Gauss_mesh(Nq_split)

#------------------------------------------------
#define function to give back state of calculation
@njit
def progressbar(i_in,imax,state):
    '''
    state=np.array([False,False,...])
    '''
    r=int(round(i_in*100/imax))
    for i in range(state.size):
        if state[i]==False and r>=round((i+1)*100/state.size):
            print(str(i+1)+'/'+str(state.size))
            state[i]=True
    return state 
#------------------------------------------------
#define function for onshell momentum
@njit
def q_0_func(E):
    q0 = np.array([np.sqrt(2*mu[0]*(E-Delta[0])+0j),np.sqrt(2*mu[1]*(E-Delta[1])+0j)])
    return q0.real+0j
#define function for COMPLEX onshell momentum
@njit
def q_0_func_complex(E):
    q0 = np.array([np.sqrt(2*mu[0]*(E-Delta[0])+0j),np.sqrt(2*mu[1]*(E-Delta[1])+0j)])
    return q0 #imaginary below threshold 
#------------------------------------------------
#define function for channel dependent two-body phase space
@njit
def two_body_phase_space1(E):
    E=E-Delta[0]
    return np.sqrt(2*mu[0]**3*E)/(4*np.pi)
@njit
def two_body_phase_space2(E):
    E=E-Delta[1]
    return np.sqrt(2*mu[1]**3*E)/(4*np.pi)
#------------------------------------------------
#define function for double factorial
@njit
def doublefactorial(n):
     if n in (0, 1):
         return 1
     else:
         return n * doublefactorial(n-2)
#define functions for Gaussian basis and its normalization
@njit
def r_n(n):
    return r_1*np.power(r_n_max/r_1,(n-1)/(n_max-1))

@njit
def nu_n(n):
    return 1/(r_n(n)**2)

@njit
def N_nl(n,l):
    return np.sqrt(2**(l+2)*(2*nu_n(n))**(l+3/2)/(np.sqrt(np.pi)*doublefactorial(2*l+1)))
#------------------------------------------------
#define functions for integrand of LSE (channel dependent)
@njit
def integrand1(q,dE):
    return 1/(dE-q**2/(2*mu[0])) 
@njit
def integrand2(q,dE):
    return 1/(dE-q**2/(2*mu[1]))
#------------------------------------------------
#------------------------------------------------
#define function to return wavefuction in position space
@njit
def psi_nlm(r,n,l,m,theta,phi):
    psi = 0
    for i in range(n_max):
        n_prime = i+1
        psi += c_solution[:,n][i]*N_nl(n_prime,l)*r**(l)*np.exp(-nu_n(n_prime)*r**2)*special.sph_harm(m,l,theta,phi)
    return psi

#define function in order to return psi with correct sign
@njit
def psi_nlm_cs(r,n,l,m,theta,phi):
    if psi_nlm(0.0+1e-4,n,l,m,theta,phi)<0.0:
        return -psi_nlm(r,n,l,m,theta,phi)
    else:
        return psi_nlm(r,n,l,m,theta,phi)

#define function for normalization of wavefunction in momentum space
@njit
def N_n0_p(n_prime):
    return np.sqrt(np.sqrt(2)/(np.pi**2)*(nu_n(n_prime)/np.pi)**(3/2))

#define function for wave function in momentum space for l=m=0:
@njit
def psi_n00_p(p,n):
    psi_p = np.zeros((n_max),dtype=np.complex128)
    for i in range(n_max):
        n_prime = i+1
        psi_p[i]=c_solution[:,n][i]*N_n0_p(n_prime)*(np.pi/nu_n(n_prime))**(3/2)*np.exp(-p**2/(4*nu_n(n_prime)))*special.sph_harm(0,0,0,0)
    
    if psi_nlm(0.0+1e-4,n,0,0,0,0)<0.0:
        return -np.sum(psi_p)
    else:
        return np.sum(psi_p)

#define function for normalization of wavefunction in momentum space for l=1,m=0
@njit
def N_n1_p(n_prime):
    return np.sqrt(4*np.sqrt(2)/(9*np.pi**(7/2))*(nu_n(n_prime))**(5/2))

#define function for wave function in momentum space for l=1,m=0:
@njit
def psi_n10_p(p,n):
    psi_p = 0.+0.j
    #p=np.abs(p)
    for i in range(n_max):
        n_prime = i+1
        psi_p += c_solution[:,n][i]*N_n1_p(n_prime)*np.sqrt(3)*np.pi*1j*1/(4*nu_n(n_prime)**(5/2))*p*np.exp(-p**2/(4*nu_n(n_prime)))
    #ensure correct sign
#    if psi_nlm(0.0+1e-4,n,1,0,0,0)<0.0:
#        return -psi_p
#    else:
#        return psi_p
    return psi_p

#define function for 2d density plot
def simple_plot2d(x,y,z,s_th=False,s_Epoles=False,v_max=2,v_min=-2,vlabel="",xl="Re $\sqrt{s}$ [GeV]",yl="Im $\sqrt{s}$ [GeV]"):
    

    fig, ax = plt.subplots()
    
    c=ax.pcolormesh(x, y, np.transpose(z), cmap='bwr', vmin=v_min, vmax=v_max,shading='auto')

    ax.set(xlabel=xl)
    ax.set(ylabel=yl)

    ax.axis([x.min(), x.max(), y.min(), y.max()])
    
    cbar=fig.colorbar(c, ax=ax)
    cbar.set_label(vlabel, rotation=90)

    if np.any(s_th)!=False:
        for i in range(s_th.size):
            ax.axvline(x=np.real(s_th[i]),linestyle="dotted",c="red")

    if np.any(s_Epoles)!=False:
        for i in range(s_Epoles.size):
            ax.axvline(x=np.real(s_Epoles[i]),linestyle="dotted",c="black")

    plt.tight_layout()
    #plt.savefig('Pole_'+string_pole+'.pdf')  
    ax.grid(True)
    plt.show()

#define sign function which returns +1 for 0 input
@njit
def sign0(E):
    if E!=0.0:
        return np.sign(E)
    else:
        return 1.0