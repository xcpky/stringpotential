#ifndef CONSTANTS_H
#define CONSTANTS_H
#define f_pi (0.092) 
#define g_pi (0.5704) 

#define m_pi (0.138039407) 
#define m_K (0.498) 
#define m_eta (0.548) 

#define m_B (5.27934) 
#define m_B_star (5.32471) 
#define m_B_s (5.36692) 
#define m_B_star_s (5.4154) 

#define m11 (m_B)
#define m12 (m_B_star)

#define m21 (m_B_s)
#define m22 (m_B_star_s)

#define m_Xb13P (10.5134-(m11+m12))
#define m_Xb12P (10.25546-(m11+m12))
#define m_Xb11P (9.89278-(m11+m12))

#define a_l (0.06426*5.068) 
#define g (0.898794378386677/a_l)
#define g1 (g*0.014926616931653945)
#define g2 (g*0.006467550544943349)
#define delta0 0
#define delta1 (m21 + m22 - m11 - m12)
#define mu0 (m11 * m12 / (m11 + m12))
#define mu1 (m21 * m22 / (m21 + m22))
#define m_c 1.85
#define a_Cornell 1.95


#endif // !DEBUG
