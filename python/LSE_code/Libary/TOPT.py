import numpy as np

from constants import *


#----------------------------------------------------------------------
#use Time Ordered Perturbation Theory (TOPT) to calculate OME potential
#----------------------------------------------------------------------
def omega(p,k,i,j):
    if i==1 and j==1:
        return 2*m11+(p**2+k**2)/(2*m11)
    if i==1 and j==2:
        return (m11+k**2/(2*m11))+(m21+p**2/(2*m21))
    if i==2 and j==1:
        return (m21+k**2/(2*m21))+(m11+p**2/(2*m11))
    if i==2 and j==2:
        return 2*m21+(p**2+k**2)/(2*m21)
    else:
        print('Wrong channel indices!')

def omega_prime(p,k,i,j):
    if i==1 and j==1:
        return 2*m12+(p**2+k**2)/(2*m12)
    if i==1 and j==2:
        return (m12+k**2/(2*m12))+(m22+p**2/(2*m22))
    if i==2 and j==1:
        return (m22+k**2/(2*m22))+(m12+p**2/(2*m12))
    if i==2 and j==2:
        return 2*m22+(p**2+k**2)/(2*m21)
    else:
        print('Wrong channel indices!')

def curlO_TOPT(E,p,k,m,i,j):
    if abs(p)<=1e-8:
        p+=1e-6
    if abs(k)<=1e-8:
        k+=1e-6
    return 4*g_pi**2/(f_pi**2)*(-1/(4*p*k)*(np.log((E-(m+(p-k)**2/(2*m))-omega(p,k,i,j))/(E-(m+(p+k)**2/(2*m))-omega(p,k,i,j)))+np.log((E-(m+(p-k)**2/(2*m))-omega_prime(p,k,i,j))/(E-(m+(p+k)**2/(2*m))-omega_prime(p,k,i,j)))))

#def function for OME potentials in TOPT
##def function for V_{11}:
def V_11_OME_TOPT(E,p,k):
    V_11 = (-3)*(3*curlO_TOPT(E,p,k,m_pi,1,1)+1/3*curlO_TOPT(E,p,k,m_eta,1,1)) #-3 factor to enforce I=0 - comes only in non strange channel
    return V_11
##def function for V_{21}:
def V_21_OME_TOPT(E,p,k):
    V_21 = 2**(3/2)*curlO_TOPT(E,p,k,m_K,2,1)
    return V_21
##def function for V_{12}:
def V_12_OME_TOPT(E,p,k):
    V_12 = 2**(3/2)*curlO_TOPT(E,p,k,m_K,1,2)
    return V_12
##def function for V_{22}:
def V_22_OME_TOPT(E,p,k):
    V_22 = 2/3*curlO_TOPT(E,p,k,m_eta,2,2)
    return V_22

def oV(E):
    dE = E - Delta + 0j
    q0 = np.sqrt(2*mu[0]*dE[0])
    return V_11_OME_TOPT(E, q0, q0)
#----------------------------------------------------------------------
#define functions for contact term
#----------------------------------------------------------------------
def ContactTerm_11(g_C):
    return (-3)*2*g_C

def ContactTerm_12(g_C):
    return (1/2)**(3/2)*4*g_C

def ContactTerm_22(g_C):
    return g_C
