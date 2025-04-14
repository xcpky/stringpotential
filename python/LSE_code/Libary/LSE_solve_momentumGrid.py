import numpy as np

from Libary.constants import *
from Libary.assistantFunc import *
from Libary.TOPT import *
from Libary.Tower import *

xq,wq=Gauss_Leg_mapping(qGauss,wGauss,0.,Lam)

@njit
def Vmat(l_E):
    E=E_list[l_E]

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq
    qval1[Nq]=q_0_func(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq
    qval2[Nq]=q_0_func(E)[1]

    vmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    for i in range(Nq+1):
        for j in range(Nq+1):
            #vmat[i,j]=Tower_psiFunc_11(qval1[i],qval1[j],l_E)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=Tower_psiFunc_12(qval2[i],qval1[j],l_E)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=Tower_psiFunc_12(qval1[i],qval2[j],l_E)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=Tower_psiFunc_22(qval2[i],qval2[j],l_E)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            vmat[i,j]=Tower_11(i,j,l_E)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            vmat[i+nn1,j]=Tower_12(j,i,l_E)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            vmat[i,j+nn1]=Tower_12(i,j,l_E)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            vmat[i+nn1,j+nn1]=Tower_22(i,j,l_E)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

    return vmat

@njit
def Vmat_psiFunc(l_E):
    E=E_list[l_E]

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq
    qval1[Nq]=q_0_func(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq
    qval2[Nq]=q_0_func(E)[1]

    vmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    for i in range(Nq+1):
        for j in range(Nq+1):
            vmat[i,j]=Tower_psiFunc_11(qval1[i],qval1[j],l_E)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            vmat[i+nn1,j]=Tower_psiFunc_12(qval2[i],qval1[j],l_E)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            vmat[i,j+nn1]=Tower_psiFunc_12(qval1[i],qval2[j],l_E)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            vmat[i+nn1,j+nn1]=Tower_psiFunc_22(qval2[i],qval2[j],l_E)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

    return vmat

@njit
def Vmat_ij(i,j,l_E):
    E=E_list[l_E]

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq
    qval1[Nq]=q_0_func(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq
    qval2[Nq]=q_0_func(E)[1]

    vmat=np.zeros((2,2),dtype=np.complex128)
    #nn1=Nq+1

    vmat[0,0]=Tower_11(i,j,l_E)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
    vmat[1,0]=Tower_12(j,i,l_E)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
    vmat[0,1]=Tower_12(i,j,l_E)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
    vmat[1,1]=Tower_22(i,j,l_E)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

    return vmat

@njit
def Vmat_psiFunc_ij(i,j,l_E):
    E=E_list[l_E]

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq
    qval1[Nq]=q_0_func(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq
    qval2[Nq]=q_0_func(E)[1]

    vmat=np.zeros((2,2),dtype=np.complex128)
    #nn1=Nq+1

    vmat[0,0]=Tower_psiFunc_11(qval1[i],qval1[j],l_E)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
    vmat[1,0]=Tower_psiFunc_12(qval2[i],qval1[j],l_E)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
    vmat[0,1]=Tower_psiFunc_12(qval1[i],qval2[j],l_E)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
    vmat[1,1]=Tower_psiFunc_22(qval2[i],qval2[j],l_E)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

    return vmat

@njit
def Gmat(l_E):
    E=E_list[l_E]
    dE1=E-Delta[0]
    dE2=E-Delta[1]

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq
    qval1[Nq]=q_0_func(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq
    qval2[Nq]=q_0_func(E)[1]

    gmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    sub1=0.+0j
    sub2=0.+0j
    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq[i]*integrand1(qval1[i],dE1)
        sub1+=wq[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq[i]*integrand2(qval2[i],dE2)
        sub2+=wq[i]*integrand2(qval2[i],dE2)
    
        if dE1>=0.0:
            gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]*(mu[0]*(-1j*np.pi+np.log(np.abs((Lam+qval1[Nq])/(Lam-qval1[Nq]))))-qval1[Nq]*sub1)
        else:
            gmat[Nq,Nq]=0.0
        if dE2>=0.0:
            gmat[Nq+nn1,Nq+nn1]=1/(2*np.pi**2)*qval2[Nq]*(mu[1]*(-1j*np.pi+np.log(np.abs((Lam+qval2[Nq])/(Lam-qval2[Nq]))))-qval2[Nq]*sub2)
        else:
            gmat[Nq+nn1,Nq+nn1]=0.0

    return gmat

@njit
def Tmat(l_E):
    vmat=Vmat(l_E)
    gmat=Gmat(l_E)

    tmat=np.linalg.inv(np.identity(2*Nq+2)-vmat@gmat)@vmat

    return tmat

@njit
def LSEpoles(l_E):
    vmat=Vmat(l_E)
    gmat=Gmat(l_E)

    return np.linalg.det(np.identity(2*Nq+2)-vmat@gmat)

@njit
def Tmat_onshell(T_grid):
    out=np.zeros((2,2,T_grid[0,0,:].size),dtype=np.complex128)
    
    nn1=Nq+1
    T11=T_grid[Nq,Nq,:]
    T12=T_grid[Nq,Nq+nn1,:]
    T21=T_grid[Nq+nn1,Nq,:]
    T22=T_grid[Nq+nn1,Nq+nn1,:]
    
    out[0,0,:]=T11
    out[0,1,:]=T12
    out[1,0,:]=T21
    out[1,1,:]=T22
    
    return out