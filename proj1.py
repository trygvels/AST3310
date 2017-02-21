#r[1/KGs]
#Q[J]
#EPS[J/KGs]

from numpy import *

#-------------------------Constants------------------------------
rho     = 1.62e5      #[kg/m^3] Solar core
T9      = 1.57e7*1e-9   #[10^9K] Solar core
N_A     = 6.0221409e23    #Dim.less
m_u     = 1.6605e-27       # [J]
m2j     = 1.6022e-13      # Conversion factor from MeV to J

#Mass fractions 
X    = 0.7
Y_3  = 1e-10
Y    = 0.29
Y_4  = Y-Y_3
Z    = 0.01
Z_Li = 10e-13
Z_Be = 10e-13
fractions = array([0, X, Y_3, Y_4, Z_Li, Z_Be]) # e, p, 3He, 4He, 7Li, 7Be

# Converting reaction energies to Joules
Q_pp = m2j*(0.15 + 5.49 + 1.02) # combined first two steps
Q_33 = m2j*12.86
Q_34 = m2j*1.59
Q_e7Be = m2j*0.05
Q_p7Li = m2j*17.35

# Number of particles for respective elements (See "fractions")
N_Particles = array([1., 1., 3., 4., 7., 7.]) #Number of particles in core. N_Particles[0]=1 to avoid division by 0

#-------------------------Lambdas-------------------------------- 
#All lambdas
lambdas = zeros((6,6))

#Lambda pp
lambda_pp = 4.01e-15 *T9**(-2./3)*exp(-3.380*T9**(-1./3))*(1+0.123*T9**(1./3)+1.09*T9**(2./3)+0.938*T9)
lambdas[1,1] = lambda_pp/(N_A*1e6) #Converting from cm3 to m3

#Lambda 33
lambda_33 = 6.04e10*T9**(-2./3)*exp(-12.276*T9**(-1./3))*(1 + 0.034*T9**(1./3) - 0.522*T9**(2./3) - 0.124*T9 + 0.353*T9**(4./3) + 0.213*T9**(-5./3))
lambdas[2,2] = lambda_33/(N_A*1e6)

#Lambda 34
T9_star = T9 / (1 + 4.95e-2 * T9)
lambda_34 = 5.61e6*T9_star**(5./6)*T9**(-3./2)*exp(-12.826*T9_star**(-1./3))
lambdas[2,3] = lambda_34/(N_A*1e6)

#Lambda e7Be
lambda_e7Be = 1.34e-10*T9**(-1./2)*(1 - 0.537*T9**(1./3) + 3.86*T9**(2./3) + 0.0027/T9*exp(2.515e-3/T9))
lambdas[0,5] = lambda_e7Be/(N_A*1e6)

#Lambda p7Li
T9_star =  T9 / (1. + 0.759*T9)
lambda_p7Li = 1.096e9*T9**(-2./3)*exp(-8.472*T9**(-1./3)) - 4.830e8 *T9_star**(5./6)*T9**(-3./2)*exp(-8.472*T9_star**(-1./3)) + 1.06e10*T9**(-3./2)*exp(-30.442/T9)
lambdas[1,4] = lambda_p7Li/(N_A*1e6)

#------------------------Calculate Number density----------------------------

n = (rho/m_u)*divide(fractions,N_Particles) # Disregarding n[0] until next line	
n[0]= n[1]+2*(n[2]+n[3])                    # n_e (Free electrons, because of ionization. No metals)

#------------------------Calculate r_ik----------------------------
r_ik=zeros((6,6)) 
for i in range(6):                  #Calculating r_ik
		for k in range(6):
			if i==k:  #Delta function
				d_ik=1
			else:
				d_ik=0
			r_ik[i,k] = lambdas[i,k]*(n[i]*n[k])/(rho*(1+d_ik)) #Lambdas makes sure this is 0, whenever we dont have a value.

#------Making Sure no step consumes more than available--------------------
if r_ik[0,5]>r_ik[2,3]:
	r_ik[0,5]=r_ik[2,3]
if r_ik[1,4]>r_ik[0,5]: 
	r_ik[1,4]=r_ik[0,5]

#------- Printing Energies? ---------------
print "r_pp (Q_pp+Q_Dp) rho = %.2e" % (r_ik[1,1]*Q_pp*rho)
print "r_33 Q_33 rho        = %.2e" % (r_ik[2,2]*Q_33*rho)
print "r_34 Q_34 rho        = %.2e" % (r_ik[2,3]*Q_34*rho)
print "r_e7Be Q_e7Be rho    = %.2e" % (r_ik[0,5]*Q_e7Be*rho)
print "r_p7Li Q_p7Li rho    = %.2e" % (r_ik[1,4]*Q_p7Li*rho)



