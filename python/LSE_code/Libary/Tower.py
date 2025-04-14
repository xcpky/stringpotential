import numpy as np

from Libary.constants import *

from Libary.assistantFunc import *

from Libary.data_saver import *

psi_mat_q1_ana = loader_3d('Grid/','psi_mat_q1_Nq'+str(Nq)+'E_listSize_'+str(E_list_length))
psi_mat_q2_ana = loader_3d('Grid/','psi_mat_q2_Nq'+str(Nq)+'E_listSize_'+str(E_list_length))

@njit
def Tower_11_complex(i,j,l_E_real,l_E_imag):

    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    
    R_11=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_11[k]=-g1**2/(E-E_solution[k])*psi_mat_q1_ana[i,l_E_real,k]*np.conj(psi_mat_q1_ana[j,l_E_real,k])
    
    return np.sum(R_11)+0j


@njit
def Tower_12_complex(i,j,l_E_real,l_E_imag):
    
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    
    R_12=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_12[k]=-g1*g2/(E-E_solution[k])*psi_mat_q1_ana[i,l_E_real,k]*np.conj(psi_mat_q2_ana[j,l_E_real,k])
    
    return np.sum(R_12)+0j

@njit
def Tower_22_complex(i,j,l_E_real,l_E_imag):
    
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    R_22=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_22[k]=-g2**2/(E-E_solution[k])*psi_mat_q2_ana[i,l_E_real,k]*np.conj(psi_mat_q2_ana[j,l_E_real,k])

    return np.sum(R_22)+0j

#--------------------------------------------------------
#only real energy values
#--------------------------------------------------------

@njit
def Tower_11(i,j,l_E):

    E=E_list[l_E]
    
    R_11=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_11[k]=-g1**2/(E-E_solution[k])*psi_mat_q1_ana[i,l_E,k]*np.conj(psi_mat_q1_ana[j,l_E,k])
    
    return np.sum(R_11)+0j


@njit
def Tower_12(i,j,l_E):

    E=E_list[l_E]
    
    R_12=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_12[k]=-g1*g2/(E-E_solution[k])*psi_mat_q1_ana[i,l_E,k]*np.conj(psi_mat_q2_ana[j,l_E,k])
    
    return np.sum(R_12)+0j

@njit
def Tower_22(i,j,l_E):

    E=E_list[l_E]

    R_22=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_22[k]=-g2**2/(E-E_solution[k])*psi_mat_q2_ana[i,l_E,k]*np.conj(psi_mat_q2_ana[j,l_E,k])

    return np.sum(R_22)+0j

#--------------------------------------------------------
#evaluate psi directly
#--------------------------------------------------------
@njit
def Tower_psiFunc_11(p,pprime,l_E):

    E=E_list[l_E]
    
    R_11=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_11[k]=-g1**2/(E-E_solution[k])*psi_n10_p(p,k)*np.conj(psi_n10_p(pprime,k))

    return np.sum(R_11)+0j


@njit
def Tower_psiFunc_12(p,pprime,l_E):

    E=E_list[l_E]

    R_12=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_12[k]=-g1*g2/(E-E_solution[k])*psi_n10_p(p,k)*np.conj(psi_n10_p(pprime,k))

    return np.sum(R_12)+0j

@njit
def Tower_psiFunc_22(p,pprime,l_E):

    E=E_list[l_E]

    R_22=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_22[k]=-g2**2/(E-E_solution[k])*psi_n10_p(p,k)*np.conj(psi_n10_p(pprime,k))

    return np.sum(R_22)+0j

#---------------------------------------------------
#define tower without wavefunctions for testing
#---------------------------------------------------

@njit
def Tower_11_complex_test(l_E_real,l_E_imag):

    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    
    R_11=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_11[k]=-g1**2/(E-E_solution[k])
    
    return np.sum(R_11)+0j


@njit
def Tower_12_complex_test(l_E_real,l_E_imag):
    
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]
    
    R_12=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_12[k]=-g1*g2/(E-E_solution[k])
    
    return np.sum(R_12)+0j

@njit
def Tower_22_complex_test(l_E_real,l_E_imag):
    
    E=E_grid_real[l_E_real]+1j*E_grid_imag[l_E_imag]

    R_22=np.zeros((NtowerUsed),dtype=np.complex128)
    
    for k in range(NtowerUsed):
        R_22[k]=-g2**2/(E-E_solution[k])

    return np.sum(R_22)+0j