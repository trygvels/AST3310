from numpy import array, loadtxt, delete,log10, pi, zeros, exp, divide
from matplotlib.pyplot import plot, show, subplot, title, xlabel,ylabel
#=============================================CONSTANTS============================================================

N_A     = 6.0221409e23    # Dim.less
m_u     = 1.6605e-27      # [kg]
m2j     = 1.6022e-13      # Conversion factor from MeV to J
k_b     =  1.38064852e-23 #[m^2 kg/(s^2 K)]
sigma   = 5.670367e-8     #[W/(m^2 K^4)]
c       = 299792458       #[m/s] 
a       = 4*sigma/c       #[kg/(m^2 s^1)]
G       = 6.674e-11       #[Nm^2/kg^2]

L_      = 3.847e26        #[W]  solar
R_      = 6.96e8          #[m]  solar
M_      = 1.989e30        #[kg] solar
rho_avg = 1.408e3         #[kg/m^3] average density of the sun
rho_0   = 5.1*rho_avg     #[kg/m^3] Initial density
T_0     = 5.7e6           #[K] Initial Temp
L_0     = 1*L_            #[J/s] Initial Luminosity 
R_0     = 0.72*R_         #[m]Total Radius 
M_0     = 1*M_            #[kg] Total mass
P_0     = 3.13e14          #Pressure

#Mass fractions 
X    = 0.7
Y_3  = 1e-10 #He3
Y = 0.29
Y_4  = Y #He4
Z    = 0.01
Z_Li = 10e-13
Z_Be = 10e-13
fractions = array([0, X, Y_3, Y_4, Z_Li, Z_Be]) # e, p, 3He, 4He, 7Li, 7Be

# mu_0 = 1/sum(Particles provided pr nucleus x massfraction/nuclear particles)
mu_0 = 1./(2*X+Y_3+3/4*Y+2/7*Z_Li)
#mu_0 = 1/(X + Y/2 + Y_3/3 + Y_4/4 + (Z_Be + Z_Li)/7)

# Converting reaction energies to Joules
Q_pp = m2j*(0.15 + 5.49 + 1.02) # combined first two steps
Q_33 = m2j*12.86
Q_34 = m2j*1.59
Q_e7Be = m2j*0.05
Q_p7Li = m2j*17.35

# Number of particles for respective elements (See "fractions")
N_Particles = array([1., 1., 3., 4., 7., 7.]) #Number of particles in core. N_Particles[0]=1 to avoid division by 0


#=============================================ENERGY PRODUCTION CALCULATION============================================================
def eps(rho,T):
	T9      = T*1e-9  
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
				if i==k:                #Delta function
					d_ik=1
				else:
					d_ik=0
				r_ik[i,k] = lambdas[i,k]*(n[i]*n[k])/(rho*(1+d_ik)) #Lambdas makes sure this is 0, whenever we dont have a value.

#------Making Sure no step consumes more than available--------------------
	if r_ik[0,5]>r_ik[2,3]:
		r_ik[0,5]=r_ik[2,3]
	if r_ik[1,4]>r_ik[0,5]: 
		r_ik[1,4]=r_ik[0,5]
#------------------Calculating Epsilon------------------------------------------
	eps = (r_ik[1,1]*Q_pp  +  r_ik[2,2]*Q_33  +  r_ik[2,3]*Q_34  +  r_ik[0,5]*Q_e7Be  +  r_ik[1,4]*Q_p7Li)
	return eps
#=============================================READING OPACITY.TXT============================================================
#Line 1 log(R) R=rho/(T/10^-6)
#Column 1 log(T) 
#rest: log(k) 

data= loadtxt('opacity.txt', skiprows=1) #reading text file, skipping first row

logR = open('opacity.txt','r').readline()   #reading first line (logR)
logR = logR.split()[1:]                     #Splitting elements in logR
logR = array([float(i) for i in logR])   #Converting elements to float

logT = data[:,0]                            #Saving first column (logT)
logK= delete(data, 0,axis=1)             #Rest of values are logK (deleting column 1)

#=============================================CALCULATING KAPPA============================================================


def findI(x,logx):      #Finding closest indices
    if logx[-1]>=x>=logx[0]:
       ic =  min(range(len(logx)), key=lambda i: abs(logx[i]-x))
       if x>=logx[ic]:  #Find two closest indices
            i=ic        #T is between i and j
            j=ic+1
       else:
            i=ic-1
            j=ic
       return i,j
    else:
        print "Value outside of boundry: Extrapolating."
        if x<logx[0]:  #Find two closest indices
            i=0        #T is between i and j
            j=1
        else:
            i=-2
            j=-1
    return i,j 
             
def kappa(rho, T):
	rho = rho*10**-3
	R = rho/(T*1e-6)**3
	R = log10(R)   
	T = log10(T)
	T_i,T_j=findI(T,logT)
	R_i,R_j=findI(R,logR)
	#Interpolation algorithm
	ki = logK[T_i][R_i]+(logK[T_j][R_i]-logK[T_i][R_i])/(logT[T_j]-logT[T_i])*(T-logT[T_i]) #k(T,r_i)
	kj = logK[T_i][R_j]+(logK[T_j][R_j]-logK[T_i][R_j])/(logT[T_j]-logT[T_i])*(T-logT[T_i]) #k(T,r_j)
	k  = ki+(kj-ki)/(logR[R_j]-logR[R_i])*(R-logR[R_i])
	kappa = 10**(k*1e-1) #Convert to SI and remove log
	return kappa
#=============================================DEFINING EQ. OF STATE ============================================================
def rho_(P,T):                                 #Equation of state Rho
	rho = (P - a/3*T**4)*mu_0*m_u/(k_b*T)
	return rho

def P_(rho,T):                                 #Equation of state Pressure
	P =  rho/(mu_0*m_u)*k_b*T+a/3*T**4
	print mu_0, m_u, k_b, a
	return P


#=============================================SOLVING DIF. EQ.============================================================

N = 10001 #Number of steps
dm = 1e26   #Step size

#-----------------------Creating empty arrays and defining initial conditions------------------------------

print rho_(P_0,T_0)
print P_(rho_0,T_0)

r = zeros(N)  
r[0] = R_0

P = zeros(N)
P[0]=P_0 #P_(rho_0,T_0)

L = zeros(N)
L[0]=L_0

rho = zeros(N)
rho[0]=rho_0

T = zeros(N)
T[0]=T_0

m = zeros(N)
m[0]=M_0

#------------------------Calculation Loop--------------------------------------------
counter = 0
for i in range(N-1): #Using the forward euler method to solfe dif. Eq.
	m[i+1] = m[i] - dm
	r[i+1] = r[i] - 1./(4*pi*r[i]**2*rho[i])*dm
	P[i+1] = P[i] + G*m[i]/(4*pi*r[i]**4)*dm
	L[i+1] = L[i] - eps(rho[i],T[i])*dm
	T[i+1] = T[i] + 3.*kappa(rho[i],T[i])*L[i]/(256.*pi**2.*sigma*r[i]**4.*T[i]**3.)*dm
	rho[i+1]= rho_(P[i+1],T[i+1])
	if m[i] < 0 or r[i] < 0 or L[i] < 0 or T[i] < 0 or rho[i] < 0: #Breaking if values go to zero
		break
	if counter == 50:          #Printing values for every 5th iteration
		print 'M: %.3g, R: %.3g, P: %.3g, L: %.3g, T: %.3g, D: %.3g' % (m[i],r[i],P[i],L[i],T[i],rho[i])		
		counter = 0
	counter +=1

#=============================================PLOTTING============================================================
subplot(2,2,1)
title('Radius vs. mass')
xlabel('M/M_sun')
ylabel('R/R_sun')
plot(m/M_,r/R_)

subplot(2,2,2)
title('Luminosity vs. mass')
xlabel('M/M_sun')
ylabel('L/L_sun')
plot(m/M_,L/L_)

subplot(2,2,3)
title('Temperature vs. mass')
xlabel('M/M_sun')
ylabel('T[K]')
plot(m/M_,T)

subplot(2,2,4)
title('Density vs. mass')
xlabel('M/M_sun')
ylabel('rho/rho_sun_avg')
plot(m/M_,rho/rho_avg)

show()
    
    
