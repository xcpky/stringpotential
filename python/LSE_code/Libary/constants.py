import numpy as np
from numba import njit

#def constants

#pion decay constant
f_pi = 0.092 #GeV
#and coupling
g_pi = 0.5704 #*2.0 #GeV

#Goldstone boson masses
m_pi = 0.138039407 #GeV
m_K = 0.498 #GeV
m_eta = 0.548 #GeV

#heavy meson masses
m_B = 5.27934 #GeV
m_B_star = 5.32471 #GeV
m_B_s = 5.36692 #GeV
m_B_star_s = 5.4154 #GeV

#asign masses to channels
m11 = m_B
m12 = m_B_star

m21 = m_B_s
m22 = m_B_star_s

#masses of bound states in 1++ channel
m_Xb13P=10.5134-(m11+m12)
m_Xb12P=10.25546-(m11+m12)
m_Xb11P=9.89278-(m11+m12)

#use couplings from D200 fit
#Lattice spacing
a_l = 0.06426*5.068 #GeV^(-1)
#pion mass dependet scale factor
g = 0.898794378386677/a_l#*1e1 #increase coupling just for illustration
#channel dependent couplings
g1 = g*0.014926616931653945
g2 = g*0.006467550544943349

#define reduced mass
mu = np.array([m11*m12/(m11+m12),m21*m22/(m21+m22)])
#and energy shift for channel 1 and 2
Delta = np.array([m11+m12,m21+m22])-(m11+m12) #m_1^0=m11, m_2^0=m12

#define constants for evaluation of wave functions (Rayleigh-Ritz principle)
n_max = 60
#r_1 = 0.02 #fm
#r_n_max = 30 #fm
r_1 = 0.02*5.068 #GeV^-1
r_n_max = 60*5.068 #GeV^-1

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#constants for solving the LSE
#cutoff
Lam=4.0 #GeV
#number of Gauss points
Nq = 40
#number of states prepared in file
NtowerFile = 30
#number of states used
NtowerUsed = 1
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#introduce width of energy poles of tower
E_width=0.01j #GeV
#load solution of generalized eigenvalue problem
# E_solution = np.loadtxt('E_solution_GaussBasis_Cornell_chib1_fit.txt')+E_width
# c_solution = np.loadtxt('c_solution_GaussBasis_Cornell_chib1_fit.txt')
#E_solution = np.loadtxt('E_solution_GaussBasis_Cornell_30_V0_04.txt')
#c_solution = np.loadtxt('c_solution_GaussBasis_Cornell_30_V0_04.txt')

#shifting zeros eigenvalue when investigating the seperable case
#E_solution[0]=0.4 #GeV

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#define energy interval for which we will scan Det and T
#E_upperBound=E_solution[0]+0.2
#E_lowerBound=E_solution[0]-0.2
E_upperBound=Delta[1]+0.5
E_lowerBound=Delta[0]-1.0
E_list_length=400
E_list = np.linspace(E_lowerBound,E_upperBound,E_list_length) #GeV
#--------------------------------------------------------------------
#define energy grid for which we will scan Det and T
E_upperBound_real=Delta[1]+0.5
E_lowerBound_real=Delta[0]-1.0

E_upperBound_imag=0.2
E_lowerBound_imag=-0.2

E_grid_size_real=400
E_grid_size_imag=200
E_grid_real = np.linspace(E_lowerBound_real,E_upperBound_real,E_grid_size_real) #GeV
E_grid_imag = np.linspace(E_lowerBound_imag,E_upperBound_imag-np.abs(E_upperBound_imag-E_lowerBound_imag)/(E_grid_size_imag),E_grid_size_imag) #GeV
#--------------------------------------------------------------------
#contact term
g_C=-0.25
#constant potential for test
V_const_test=1.0
