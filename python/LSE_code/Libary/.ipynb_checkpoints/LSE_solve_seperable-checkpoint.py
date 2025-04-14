from Libary.Tower import *
from Libary.assistantFunc import *

#check correction of code via analytical solution of T matrix of the tower potential with just one state (which is therefore seperable)

#define function for G matrix
@njit
def Gmat1_test(l_E):
    E=E_list[l_E]
    dE1=E-Delta[0]+0j
    q01=q_0_func(E)[0]

    qval_grid,w_q=Gauss_Leg_mapping(qGauss,wGauss,0.,Lam)

    res=0.+0j
    sub=0.+0j
    for i in range(qval_grid.size):
        #res+=1/(2*np.pi**2)*qval_grid[i]**2*w_q[i]*np.abs(psi_mat_q1_ana[0+NtowerFile*l_E,i])**2*integrand1(qval_grid[i],dE1)
        res+=1/(2*np.pi**2)*qval_grid[i]**2*w_q[i]*np.abs(psi_mat_q1_ana[i,l_E,0])**2*integrand1(qval_grid[i],dE1)
        sub+=w_q[i]*integrand1(qval_grid[i],dE1)
    if E>=Delta[0]:
        #res+=1/(2*np.pi**2)*np.abs(psi_mat_q1_ana[0+NtowerFile*l_E,Nq])**2*q01*(mu[0]*(-1j*np.pi + np.log(np.abs(Lam+q01)/np.abs(Lam-q01)))-q01*sub)
        res+=1/(2*np.pi**2)*np.abs(psi_mat_q1_ana[Nq,l_E,0])**2*q01*(mu[0]*(-1j*np.pi + np.log(np.abs(Lam+q01)/np.abs(Lam-q01)))-q01*sub)
    
    return res

@njit
def Gmat2_test(l_E):
    E=E_list[l_E]
    dE2=E-Delta[1]+0j
    q02=q_0_func(E)[1]

    qval_grid,w_q=Gauss_Leg_mapping(qGauss,wGauss,0.,Lam)

    res=0.+0j
    sub=0.+0j
    for i in range(qval_grid.size):
        #res+=1/(2*np.pi**2)*qval_grid[i]**2*w_q[i]*np.abs(psi_mat_q2_ana[0+NtowerFile*l_E,i])**2*integrand2(qval_grid[i],dE2)
        res+=1/(2*np.pi**2)*qval_grid[i]**2*w_q[i]*np.abs(psi_mat_q2_ana[i,l_E,0])**2*integrand2(qval_grid[i],dE2)
        sub+=w_q[i]*integrand2(qval_grid[i],dE2)
    if E>=Delta[1]:
        #res+=1/(2*np.pi**2)*np.abs(psi_mat_q2_ana[0+NtowerFile*l_E,Nq])**2*q02*(mu[1]*(-1j*np.pi + np.log(np.abs(Lam+q02)/np.abs(Lam-q02)))-q02*sub)
        res+=1/(2*np.pi**2)*np.abs(psi_mat_q2_ana[Nq,l_E,0])**2*q02*(mu[1]*(-1j*np.pi + np.log(np.abs(Lam+q02)/np.abs(Lam-q02)))-q02*sub)
    
    return res

@njit
def onshellTmat_test(l_E):
    Tmat=np.zeros((2,2),dtype=np.complex128)

    G1=Gmat1_test(l_E)
    G2=Gmat2_test(l_E)

    q01=q_0_func(E_list[l_E])[0]
    q02=q_0_func(E_list[l_E])[1]

    E=E_list[l_E]

    Tmat[0,0]=g1**2*np.abs(psi_mat_q1_ana[Nq,l_E,0])**2/(E-E_solution[0]+g1**2*G1+g2**2*G2)
    Tmat[0,1]=g1*g2*psi_mat_q1_ana[Nq,l_E,0]*np.conj(psi_mat_q2_ana[Nq,l_E,0])/(E-E_solution[0]+g1**2*G1+g2**2*G2)
    Tmat[1,0]=g2*g1*psi_mat_q2_ana[Nq,l_E,0]*np.conj(psi_mat_q1_ana[Nq,l_E,0])/(E-E_solution[0]+g1**2*G1+g2**2*G2)
    Tmat[1,1]=g2**2*np.abs(psi_mat_q2_ana[Nq,l_E,0])**2/(E-E_solution[0]+g1**2*G1+g2**2*G2)

    #Tmat[0,0]=g1**2*np.abs(psi_n10_p(q01,0))**2/(E-E_solution[0]+g1**2*G1+g2**2*G2)
    #Tmat[0,1]=g1*g2*psi_n10_p(q01,0)*np.conj(psi_n10_p(q02,0))/(E-E_solution[0]+g1**2*G1+g2**2*G2)
    #Tmat[1,0]=g2*g1*psi_n10_p(q02,0)*np.conj(psi_n10_p(q01,0))/(E-E_solution[0]+g1**2*G1+g2**2*G2)
    #Tmat[1,1]=g2**2*np.abs(psi_n10_p(q02,0))**2/(E-E_solution[0]+g1**2*G1+g2**2*G2)

    return Tmat

@njit
def Det1mVG_test(l_E):
    G1=Gmat1_test(l_E)
    G2=Gmat2_test(l_E)

    E=E_list[l_E]

    return 1+g1**2/(E-E_solution[0])*G1+g2**2/(E-E_solution[0])*G2
