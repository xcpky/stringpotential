import numpy as np

from Libary.constants import *
from Libary.assistantFunc import *
from Libary.TOPT import *
from Libary.Tower import *

#xq,wq=Gauss_Leg_mapping(qGauss,wGauss,0.,Lam)

#split integration into two parts, i.e. 1st threshold to 2nd threshold, 2nd threshold to Lambda
qGaussFullInt,wGaussFullInt=Gauss_mesh(Nq)
xq,wq=Gauss_Leg_mapping(qGaussFullInt,wGaussFullInt,0.,Lam)

@njit
def Vmat_below(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    xq1=np.zeros((Nq),dtype=np.complex128)
    wq1=np.zeros((Nq),dtype=np.complex128)
    xq1[:Nq_split],wq1[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[0])
    xq1[Nq_split:],wq1[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[0],Lam)

    xq2=np.zeros((Nq),dtype=np.complex128)
    wq2=np.zeros((Nq),dtype=np.complex128)
    xq2[:Nq_split],wq2[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[1])
    xq2[Nq_split:],wq2[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[1],Lam)

    qval1=np.zeros((Nq),dtype=np.complex128)
    qval1=xq1

    qval2=np.zeros((Nq),dtype=np.complex128)
    qval2=xq2

    vmat=np.zeros((2*Nq,2*Nq),dtype=np.complex128)
    nn1=Nq

    for i in range(nn1):
        for j in range(nn1):

            vmat[i,j]=Tower_11_complex(i,j,l_E_real,l_E_imag)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            vmat[i+nn1,j]=Tower_12_complex(j,i,l_E_real,l_E_imag)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            vmat[i,j+nn1]=Tower_12_complex(i,j,l_E_real,l_E_imag)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            vmat[i+nn1,j+nn1]=Tower_22_complex(i,j,l_E_real,l_E_imag)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=Tower_11_complex_test(l_E_real,l_E_imag)
            #vmat[i+nn1,j]=Tower_12_complex_test(l_E_real,l_E_imag)
            #vmat[i,j+nn1]=Tower_12_complex_test(l_E_real,l_E_imag)
            #vmat[i+nn1,j+nn1]=Tower_22_complex_test(l_E_real,l_E_imag)

    return vmat

@njit
def Vmat_between(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    xq1=np.zeros((Nq),dtype=np.complex128)
    wq1=np.zeros((Nq),dtype=np.complex128)
    xq1[:Nq_split],wq1[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[0])
    xq1[Nq_split:],wq1[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[0],Lam)

    xq2=np.zeros((Nq),dtype=np.complex128)
    wq2=np.zeros((Nq),dtype=np.complex128)
    xq2[:Nq_split],wq2[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[1])
    xq2[Nq_split:],wq2[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[1],Lam)

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq1
    qval1[Nq]=q_0_func(E_grid_real[l_E_real])[0]

    qval2=np.zeros((Nq),dtype=np.complex128)
    qval2=xq2

    vmat=np.zeros((2*Nq+1,2*Nq+1),dtype=np.complex128)
    nn1=Nq+1
    nn2=Nq

    for i in range(nn1):
        for j in range(nn1):
            vmat[i,j]=Tower_11_complex(i,j,l_E_real,l_E_imag)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i,j]=V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i,j]=Tower_11_complex_test(l_E_real,l_E_imag)
    for i in range(nn2):
        for j in range(nn1):
            vmat[i+nn1,j]=Tower_12_complex(j,i,l_E_real,l_E_imag)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j]=V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j]=Tower_12_complex_test(l_E_real,l_E_imag)
    for i in range(nn1):
        for j in range(nn2):
            vmat[i,j+nn1]=Tower_12_complex(i,j,l_E_real,l_E_imag)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=Tower_12_complex_test(l_E_real,l_E_imag)
    for i in range(nn2):
        for j in range(nn2):
            vmat[i+nn1,j+nn1]=Tower_22_complex(i,j,l_E_real,l_E_imag)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)
            #vmat[i+nn1,j+nn1]=V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)
            #vmat[i+nn1,j+nn1]=Tower_22_complex_test(l_E_real,l_E_imag)

    return vmat

@njit
def Vmat_above(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    xq1=np.zeros((Nq),dtype=np.complex128)
    wq1=np.zeros((Nq),dtype=np.complex128)
    xq1[:Nq_split],wq1[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[0])
    xq1[Nq_split:],wq1[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[0],Lam)

    xq2=np.zeros((Nq),dtype=np.complex128)
    wq2=np.zeros((Nq),dtype=np.complex128)
    xq2[:Nq_split],wq2[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[1])
    xq2[Nq_split:],wq2[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[1],Lam)

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq1
    qval1[Nq]=q_0_func(E_grid_real[l_E_real])[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq2
    qval2[Nq]=q_0_func(E_grid_real[l_E_real])[1]

    vmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    for i in range(Nq+1):
        for j in range(Nq+1):

            vmat[i,j]=Tower_11_complex(i,j,l_E_real,l_E_imag)#+V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            vmat[i+nn1,j]=Tower_12_complex(j,i,l_E_real,l_E_imag)#+V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            vmat[i,j+nn1]=Tower_12_complex(i,j,l_E_real,l_E_imag)#+V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            vmat[i+nn1,j+nn1]=Tower_22_complex(i,j,l_E_real,l_E_imag)#+V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=V_11_OME_TOPT(E,qval1[i],qval1[j])+ContactTerm_11(g_C)
            #vmat[i+nn1,j]=V_21_OME_TOPT(E,qval2[i],qval1[j])+ContactTerm_12(g_C)
            #vmat[i,j+nn1]=V_12_OME_TOPT(E,qval1[i],qval2[j])+ContactTerm_12(g_C)
            #vmat[i+nn1,j+nn1]=V_22_OME_TOPT(E,qval2[i],qval2[j])+ContactTerm_22(g_C)

            #vmat[i,j]=Tower_11_complex_test(l_E_real,l_E_imag)
            #vmat[i+nn1,j]=Tower_12_complex_test(l_E_real,l_E_imag)
            #vmat[i,j+nn1]=Tower_12_complex_test(l_E_real,l_E_imag)
            #vmat[i+nn1,j+nn1]=Tower_22_complex_test(l_E_real,l_E_imag)

    return vmat

@njit
def Gmat_below(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    dE1=E-Delta[0]
    dE2=E-Delta[1]

    xq1=np.zeros((Nq),dtype=np.complex128)
    wq1=np.zeros((Nq),dtype=np.complex128)
    xq1[:Nq_split],wq1[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[0])
    xq1[Nq_split:],wq1[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[0],Lam)

    xq2=np.zeros((Nq),dtype=np.complex128)
    wq2=np.zeros((Nq),dtype=np.complex128)
    xq2[:Nq_split],wq2[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[1])
    xq2[Nq_split:],wq2[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[1],Lam)

    qval1=np.zeros((Nq),dtype=np.complex128)
    qval1=xq1

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2=xq2

    gmat=np.zeros((2*Nq,2*Nq),dtype=np.complex128)
    nn1=Nq

    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq1[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq2[i]*integrand2(qval2[i],dE2)

    return gmat

@njit
def Gmat_between(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    dE1=E-Delta[0]
    dE2=E-Delta[1]

    xq1=np.zeros((Nq),dtype=np.complex128)
    wq1=np.zeros((Nq),dtype=np.complex128)
    xq1[:Nq_split],wq1[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[0])
    xq1[Nq_split:],wq1[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[0],Lam)

    xq2=np.zeros((Nq),dtype=np.complex128)
    wq2=np.zeros((Nq),dtype=np.complex128)
    xq2[:Nq_split],wq2[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[1])
    xq2[Nq_split:],wq2[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[1],Lam)

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq1
    qval1[Nq]=q_0_func(E_grid_real[l_E_real])[0]
    #q01c=q_0_func_complex(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2=xq2

    gmat=np.zeros((2*Nq+1,2*Nq+1),dtype=np.complex128)
    nn1=Nq+1

    sub1=0.+0j
    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq1[i]*integrand1(qval1[i],dE1)
        sub1+=wq1[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq2[i]*integrand2(qval2[i],dE2)
    
    gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]*(mu[0]*(-sign0(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval1[Nq])/(Lam-qval1[Nq]))))-qval1[Nq]*sub1)
    #gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]**2*(mu[0]/q01c*(np.log((Lam+q01c)/(Lam-q01c)))-sub1)

    return gmat

@njit
def Gmat_above(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    dE1=E-Delta[0]
    dE2=E-Delta[1]

    xq1=np.zeros((Nq),dtype=np.complex128)
    wq1=np.zeros((Nq),dtype=np.complex128)
    xq1[:Nq_split],wq1[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[0])
    xq1[Nq_split:],wq1[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[0],Lam)

    xq2=np.zeros((Nq),dtype=np.complex128)
    wq2=np.zeros((Nq),dtype=np.complex128)
    xq2[:Nq_split],wq2[:Nq_split]=Gauss_Leg_mapping(qGauss,wGauss,0.,q_0_func(E_grid_real[l_E_real])[1])
    xq2[Nq_split:],wq2[Nq_split:]=Gauss_Leg_mapping(qGauss,wGauss,q_0_func(E_grid_real[l_E_real])[1],Lam)

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq1
    qval1[Nq]=q_0_func(E_grid_real[l_E_real])[0]
    #q01c=q_0_func_complex(E)[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq2
    qval2[Nq]=q_0_func(E_grid_real[l_E_real])[1]
    #q02c=q_0_func_complex(E)[1]

    gmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    sub1=0.+0j
    sub2=0.+0j
    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq1[i]*integrand1(qval1[i],dE1)
        sub1+=wq1[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq2[i]*integrand2(qval2[i],dE2)
        sub2+=wq2[i]*integrand2(qval2[i],dE2)
    
    
    gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]*(mu[0]*(-sign0(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval1[Nq])/(Lam-qval1[Nq]))))-qval1[Nq]*sub1)
    #gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]**2*(mu[0]/q01c*(np.log((Lam+q01c)/(Lam-q01c)))-sub1)

    gmat[Nq+nn1,Nq+nn1]=1/(2*np.pi**2)*qval2[Nq]*(mu[1]*(-sign0(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval2[Nq])/(Lam-qval2[Nq]))))-qval2[Nq]*sub2)
    #gmat[Nq+nn1,Nq+nn1]=1/(2*np.pi**2)*qval2[Nq]**2*(mu[1]/q02c*(np.log((Lam+q02c)/(Lam-q02c)))-sub2)

    return gmat

@njit
def Gmat_test(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    dE1=E-Delta[0]
    dE2=E-Delta[1]

    qval1=np.zeros((Nq+1),dtype=np.complex128)
    qval1[:-1]=xq
    qval1[Nq]=q_0_func(E_grid_real[l_E_real])[0]

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq
    qval2[Nq]=q_0_func(E_grid_real[l_E_real])[1]

    gmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    sub1=0.+0j
    sub2=0.+0j
    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq[i]*integrand1(qval1[i],dE1)
        sub1+=wq[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq[i]*integrand2(qval2[i],dE2)
        sub2+=wq[i]*integrand2(qval2[i],dE2)
    
    if E.real>=Delta[0] and np.abs(E.imag)<=1e-10:
        gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]*(mu[0]*(-np.sign(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval1[Nq])/(Lam-qval1[Nq]))))-qval1[Nq]*sub1)
    else:
        gmat[Nq,Nq]=0.0+0j
    
    if E.real>=Delta[0] and np.abs(E.imag)<=1e-10:
        gmat[Nq+nn1,Nq+nn1]=1/(2*np.pi**2)*qval2[Nq]*(mu[1]*(-np.sign(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval2[Nq])/(Lam-qval2[Nq]))))-qval2[Nq]*sub2)
    else:
        gmat[Nq+nn1,Nq+nn1]=0.0+0j

    return gmat

@njit
def Gmat_complex_test(l_E_real,l_E_imag):

    gmat=np.zeros((2,2),dtype=np.complex128)

    gmat_help=Gmat_test(l_E_real,l_E_imag)
    gmat[0,0]=np.trace(gmat_help[:Nq,:Nq])
    gmat[1,1]=np.trace(gmat_help[Nq:,Nq:])

    return gmat

@njit
def Tmat_complex(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    if np.abs(E.imag)>=1e-10 or E.real<Delta[0]:
        vmat=Vmat_below(l_E_real,l_E_imag)
        gmat=Gmat_below(l_E_real,l_E_imag)
        one=np.identity(2*Nq)
    if E.real>=Delta[0] and E.real<Delta[1]:
        vmat=Vmat_between(l_E_real,l_E_imag)
        gmat=Gmat_between(l_E_real,l_E_imag)
        one=np.identity(2*Nq+1)
    else:
        vmat=Vmat_above(l_E_real,l_E_imag)
        gmat=Gmat_above(l_E_real,l_E_imag)
        one=np.identity(2*Nq+2)

    tmat=np.linalg.inv(one-vmat@gmat)@vmat

    return tmat

@njit
def LSEpoles_complex(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    if np.abs(E.imag)>=1e-10 or E.real<Delta[0]:
        vmat=Vmat_below(l_E_real,l_E_imag)
        gmat=Gmat_below(l_E_real,l_E_imag)
        one=np.identity(2*Nq)
    elif E.real>=Delta[0] and E.real<Delta[1]:
        vmat=Vmat_between(l_E_real,l_E_imag)
        gmat=Gmat_between(l_E_real,l_E_imag)
        one=np.identity(2*Nq+1)
    else:
        vmat=Vmat_above(l_E_real,l_E_imag)
        gmat=Gmat_above(l_E_real,l_E_imag)
        one=np.identity(2*Nq+2)

    return np.linalg.det(one-vmat@gmat)

@njit
def Gmat_complex(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    gmat=np.zeros((2,2),dtype=np.complex128)

    if np.abs(E.imag)>=1e-10 or E.real<Delta[0]:
        gmat_help=Gmat_below(l_E_real,l_E_imag)
        gmat[0,0]=np.trace(gmat_help[:Nq,:Nq])
        gmat[1,1]=np.trace(gmat_help[Nq:,Nq:])
    elif E.real>=Delta[0] and E.real<Delta[1]:
        gmat_help=Gmat_between(l_E_real,l_E_imag)
        gmat[0,0]=np.trace(gmat_help[:Nq+1,:Nq+1])
        gmat[1,1]=np.trace(gmat_help[Nq+1:,Nq+1:])
    else:
        gmat_help=Gmat_above(l_E_real,l_E_imag)
        gmat[0,0]=np.trace(gmat_help[:Nq+1,:Nq+1])
        gmat[1,1]=np.trace(gmat_help[Nq+1:,Nq+1:])

    return gmat

@njit
def Tmat_onshell_complex(T_grid):
    out=np.zeros((2,2,T_grid[0,0,:,0].size,T_grid[0,0,0,:].size),dtype=np.complex128)
    
    nn1=Nq+1
    T11=T_grid[Nq,Nq,:,:]
    T12=T_grid[Nq,Nq+nn1,:,:]
    T21=T_grid[Nq+nn1,Nq,:,:]
    T22=T_grid[Nq+nn1,Nq+nn1,:,:]
    
    out[0,0,:,:]=T11
    out[0,1,:,:]=T12
    out[1,0,:,:]=T21
    out[1,1,:,:]=T22
    
    return out