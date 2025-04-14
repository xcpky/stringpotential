from Libary.Tower import *
from Libary.assistantFunc import *
from Libary.LSE_solve_momentumGrid_complex import *

#check correction of code via analytical solution of T matrix of the tower potential with just one state (which is therefore seperable)

#define function for G matrix

@njit
def Gmat_wpsi_below(l_E_real,l_E_imag):
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
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq1[i]*np.abs(psi_mat_q1_ana[i,l_E_real,0])**2*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq2[i]*np.abs(psi_mat_q2_ana[i,l_E_real,0])**2*integrand2(qval2[i],dE2)

    return gmat

@njit
def Gmat_wpsi_between(l_E_real,l_E_imag):
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

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2=xq2

    gmat=np.zeros((2*Nq+1,2*Nq+1),dtype=np.complex128)
    nn1=Nq+1

    sub1=0.+0j
    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq1[i]*np.abs(psi_mat_q1_ana[i,l_E_real,0])**2*integrand1(qval1[i],dE1)
        sub1+=wq1[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq2[i]*np.abs(psi_mat_q2_ana[i,l_E_real,0])**2*integrand2(qval2[i],dE2)
    
    gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]*np.abs(psi_mat_q1_ana[Nq,l_E_real,0])**2*(mu[0]*(-np.sign(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval1[Nq])/(Lam-qval1[Nq]))))-qval1[Nq]*sub1)

    return gmat

@njit
def Gmat_wpsi_above(l_E_real,l_E_imag):
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

    qval2=np.zeros((Nq+1),dtype=np.complex128)
    qval2[:-1]=xq2
    qval2[Nq]=q_0_func(E_grid_real[l_E_real])[1]

    gmat=np.zeros((2*Nq+2,2*Nq+2),dtype=np.complex128)
    nn1=Nq+1

    sub1=0.+0j
    sub2=0.+0j
    for i in range(Nq):
        gmat[i,i]=1/(2*np.pi**2)*qval1[i]**2*wq1[i]*np.abs(psi_mat_q1_ana[i,l_E_real,0])**2*integrand1(qval1[i],dE1)
        sub1+=wq1[i]*integrand1(qval1[i],dE1)
        
        gmat[i+nn1,i+nn1]=1/(2*np.pi**2)*qval2[i]**2*wq2[i]*np.abs(psi_mat_q2_ana[i,l_E_real,0])**2*integrand2(qval2[i],dE2)
        sub2+=wq2[i]*integrand2(qval2[i],dE2)
    
    
    gmat[Nq,Nq]=1/(2*np.pi**2)*qval1[Nq]*np.abs(psi_mat_q1_ana[Nq,l_E_real,0])**2*(mu[0]*(-np.sign(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval1[Nq])/(Lam-qval1[Nq]))))-qval1[Nq]*sub1)

    gmat[Nq+nn1,Nq+nn1]=1/(2*np.pi**2)*qval2[Nq]*np.abs(psi_mat_q2_ana[Nq,l_E_real,0])**2*(mu[1]*(-np.sign(E.imag)*1j*np.pi+np.log(np.abs((Lam+qval2[Nq])/(Lam-qval2[Nq]))))-qval2[Nq]*sub2)

    return gmat

@njit
def Gmat_wpsi_complex(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    gmat=np.zeros((2,2),dtype=np.complex128)

    if np.abs(E.imag)>=1e-10 or E.real<Delta[0]:
        gmat_help=Gmat_wpsi_below(l_E_real,l_E_imag)
        gmat[0,0]=np.trace(gmat_help[:Nq,:Nq])
        gmat[1,1]=np.trace(gmat_help[Nq:,Nq:])
    elif E.real>=Delta[0] and E.real<Delta[1]:
        gmat_help=Gmat_wpsi_between(l_E_real,l_E_imag)
        gmat[0,0]=np.trace(gmat_help[:Nq+1,:Nq+1])
        gmat[1,1]=np.trace(gmat_help[Nq+1:,Nq+1:])
    else:
        gmat_help=Gmat_wpsi_above(l_E_real,l_E_imag)
        gmat[0,0]=np.trace(gmat_help[:Nq+1,:Nq+1])
        gmat[1,1]=np.trace(gmat_help[Nq+1:,Nq+1:])

    return gmat

@njit
def Det1mVG_complex(l_E_real,l_E_imag):
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    if np.abs(E.imag)>=1e-10 or E.real<Delta[0]:
        gmat_help=Gmat_wpsi_below(l_E_real,l_E_imag)
        G1=np.trace(gmat_help[:Nq,:Nq])
        G2=np.trace(gmat_help[Nq:,Nq:])
    elif E.real>=Delta[0] and E.real<Delta[1]:
        gmat_help=Gmat_wpsi_between(l_E_real,l_E_imag)
        G1=np.trace(gmat_help[:Nq+1,:Nq+1])
        G2=np.trace(gmat_help[Nq+1:,Nq+1:])
    else:
        gmat_help=Gmat_wpsi_above(l_E_real,l_E_imag)
        G1=np.trace(gmat_help[:Nq+1,:Nq+1])
        G2=np.trace(gmat_help[Nq+1:,Nq+1:])

    return 1+g1**2/(E-E_solution[0])*G1+g2**2/(E-E_solution[0])*G2
