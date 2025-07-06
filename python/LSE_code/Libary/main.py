from constants import *
from TOPT import *
import numpy as np
import matplotlib.pyplot as plt

def Gauss_mesh(n):
    '''
    return xi(n),wi(n)
    '''
    out=np.polynomial.legendre.leggauss(n)
    return out[0],out[1]

def Gauss_Leg_mapping(x,a,b):
    return (b-a)*x/2+(a+b)/2
    
def Gauss_Legendre(f,a,b,x,w,*args):
    '''
    f=function(x_int,*args)
    '''
    res=0.+0j
    x_var=Gauss_Leg_mapping(x,a,b)
    for i in range(x_var.size):
        res+=w[i]*f(x_var[i],*args)
    return (b-a)*res/2

mu = np.array([m11*m12/(m11+m12),m21*m22/(m21+m22)])

Delta = np.array([m11+m12,m21+m22])-(m11+m12) #m_1^0=m11, m_2^0=m12

def q_0_func(E):
    q0 = np.array([np.sqrt(2*mu[0]*(E-Delta[0])+0j),np.sqrt(2*mu[1]*(E-Delta[1])+0j)])
    return q0 #imaginary below threshold

#@njit
def integrand1(q,dE):
    return 1/(dE-q**2/(2*mu[0])) 
#@njit
def integrand2(q,dE):
    return 1/(dE-q**2/(2*mu[1]))

def Vmat(N_q,Lam,E,k_E,g_C):
    xi_q,wi_q=Gauss_mesh(N_q)

    qval_grid=Gauss_Leg_mapping(xi_q,0.,Lam)+0j

    vmat=np.zeros((2*N_q+2,2*N_q+2),dtype=np.complex128)
    
    q01=q_0_func(E)[0]
    q02=q_0_func(E)[1]
    
    qval1=np.zeros((N_q+1),dtype=np.complex128)
    qval1[:-1]=qval_grid
    qval1[N_q]=q01

    qval2=np.zeros((N_q+1),dtype=np.complex128)
    qval2[:-1]=qval_grid
    qval2[N_q]=q02

    nn1=N_q+1
    for i in range(N_q+1):
        for j in range(N_q+1):
            #vmat[i,j]=Tower_11_test(E)
            #vmat[i+nn1,j]=Tower_12_test(E)
            #vmat[i,j+nn1]=Tower_12_test(E)
            #vmat[i+nn1,j+nn1]=Tower_22_test(E)

            vmat[i,j]=V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            vmat[i+nn1,j]=V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            vmat[i,j+nn1]=V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            vmat[i+nn1,j+nn1]=V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=V_11_OME(qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=V_21_OME(qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=V_12_OME(qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=V_22_OME(qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=Tower_11_psi_Test(E,qval1[i],qval1[j])
            #vmat[i+nn1,j]=Tower_12_psi_Test(E,qval1[i],qval2[j])
            #vmat[i,j+nn1]=Tower_21_psi_Test(E,qval2[i],qval1[j])
            #vmat[i+nn1,j+nn1]=Tower_22_psi_Test(E,qval2[i],qval2[j])

            #vmat[i,j]=Tower_11_test(E)+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=Tower_12_test(E)+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=Tower_12_test(E)+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=Tower_22_test(E)+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=Tower_11(i,j,k_E,E)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=Tower_12(i,j,k_E,E)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=Tower_12(i,j,k_E,E)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=Tower_22(i,j,k_E,E)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)
    
    return vmat

def Gmat(N_q,Lam,E):
    xi_q,wi_q=Gauss_mesh(N_q)

    qval_grid=Gauss_Leg_mapping(xi_q,0.,Lam)+0j

    gmat=np.zeros((2*N_q+2,2*N_q+2),dtype=np.complex128)
    
    dE1=E-Delta[0]+0j
    dE2=E-Delta[1]+0j
    
    q01=q_0_func(E)[0]
    q02=q_0_func(E)[1]
    

    for i in range(N_q):
        gmat[i,i]=1/(2*np.pi**2)*qval_grid[i]**2*wi_q[i]*Lam/2./(dE1-qval_grid[i]**2/(2*mu[0]))
        gmat[i+N_q+1,i+N_q+1]=1/(2*np.pi**2)*qval_grid[i]**2*wi_q[i]*Lam/2./(dE2-qval_grid[i]**2/(2*mu[1]))
    if E>=Delta[0]:   
        gmat[N_q,N_q]=1/(2*np.pi**2)*q01*(mu[0]*(-1j*np.pi + np.log(np.abs(Lam+q01)/np.abs(Lam-q01)))-q01*Gauss_Legendre(integrand1,0.,Lam,xi_q,wi_q,dE1))
    else:
        gmat[N_q,N_q]=0.0+0j
    if E>=Delta[1]:
        gmat[2*N_q+1,2*N_q+1]=1/(2*np.pi**2)*q02*(mu[1]*(-1j*np.pi + np.log(np.abs(Lam+q02)/np.abs(Lam-q02)))-q02*Gauss_Legendre(integrand2,0.,Lam,xi_q,wi_q,dE2))
    else:
        gmat[2*N_q+1,2*N_q+1]=0.0+0j
    
    return gmat

def T_matrix(N_q,Lam,E,E_j,g_C):
    V=Vmat(N_q,Lam,E + m11 + m12,E_j,g_C)
    G=Gmat(N_q,Lam,E)
    T=np.linalg.pinv(np.identity(2*N_q+2)-V@G)@V
    #------------------------------------------
    VG_inv=np.linalg.pinv(np.identity(2*N_q+2)-V@G)
    VG=np.identity(2*N_q+2)-V@G
    if np.allclose(np.dot(VG_inv,VG),np.identity(2*N_q+2))==False:
        print('Inverting failed!','VG-1VG=',np.dot(VG_inv,VG))
    #------------------------------------------
    return T,VG,VG_inv,V,G

N_q = 40
Lam = 4
def onshellV(E):
    V=Vmat(N_q,Lam,E,0,g_C)
    return V[Nq, Nq], V[2*Nq+1, 2*Nq+1]

def onshellG(E):
    G=Gmat(N_q, Lam, E)
    return G[Nq, Nq], G[2*Nq+1, 2*Nq+1]

def onshellT(E):
    T = T_matrix(N_q, Lam, E, 0, g_C)[0]
    return T[Nq, Nq], T[2*Nq+1, 2*Nq+1]

E = np.linspace(-2.5, 0.4, 2000)
e = 0.2
V = Vmat(N_q, Lam, e, 0, g_C)
oT = [onshellT(e) for e in E]
oT11 = [oT[i][0] for i in range(len(E))]
oT22 = [oT[i][1] for i in range(len(E))]
plt.plot(E, np.abs(oT))
plt.ylim(0,1e5)
plt.savefig("ot.png",dpi=300)
# detV = [np.linalg.det(Vmat(N_q, Lam, e, 0, g_C)) for e in E]
