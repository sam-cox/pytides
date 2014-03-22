
import numpy as np
d2r, r2d = np.pi/180.0, 180.0/np.pi

#The following functions take a dictionary of astronomical values (in degrees)
#and return dimensionless scale factors for constituent amplitudes.

def f_unity(a):
	return 1.0

#Schureman equations 73, 65
def f_Mm(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = (2/3.0 - np.sin(omega)**2)*(1 - 3/2.0 * np.sin(i)**2)
	return (2/3.0 - np.sin(I)**2) / mean

#Schureman equations 74, 66
def f_Mf(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(omega)**2 * np.cos(0.5*i)**4
	return np.sin(I)**2 / mean

#Schureman equations 75, 67
def f_O1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(omega) * np.cos(0.5*omega)**2 * np.cos(0.5*i)**4
	return (np.sin(I) * np.cos(0.5*I)**2) / mean

#Schureman equations 76, 68
def f_J1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(2*omega) * (1-3/2.0 * np.sin(i)**2)
	return np.sin(2*I) / mean

#Schureman equations 77, 69
def f_OO1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(omega) * np.sin(0.5*omega)**2 * np.cos(0.5*i)**4
	return np.sin(I) * np.sin(0.5*I)**2 / mean

#Schureman equations 78, 70
def f_M2(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.cos(0.5*omega)**4 * np.cos(0.5*i)**4
	return np.cos(0.5*I)**4 / mean

#Schureman equations 227, 226, 68
#Should probably eventually include the derivations of the magic numbers (0.5023 etc).
def f_K1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	nu = d2r*a['nu'].value
	sin2Icosnu_mean = np.sin(2*omega) * (1-3/2.0 * np.sin(i)**2)
	mean = 0.5023*sin2Icosnu_mean + 0.1681
	return (0.2523*np.sin(2*I)**2 + 0.1689*np.sin(2*I)*np.cos(nu)+0.0283)**(0.5) / mean

#Schureman equations 215, 213, 204
#It can be (and has been) confirmed that the exponent for R_a reads 1/2 via Schureman Table 7
def f_L2(a):
	P = d2r*a['P'].value
	I = d2r*a['I'].value
	R_a_inv = (1 - 12*np.tan(0.5*I)**2 * np.cos(2*P)+36*np.tan(0.5*I)**4)**(0.5)
	return f_M2(a) * R_a_inv

#Schureman equations 235, 234, 71
#Again, magic numbers
def f_K2(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	nu = d2r*a['nu'].value
	sinsqIcos2nu_mean = np.sin(omega)**2 * (1-3/2.0 * np.sin(i)**2)
	mean = 0.5023*sinsqIcos2nu_mean + 0.0365
	return (0.2533*np.sin(I)**4 + 0.0367*np.sin(I)**2 *np.cos(2*nu)+0.0013)**(0.5) / mean

#Schureman equations 206, 207, 195
def f_M1(a):
	P = d2r*a['P'].value
	I = d2r*a['I'].value
	Q_a_inv = (0.25 + 1.5*np.cos(I)*np.cos(2*P)*np.cos(0.5*I)**(-0.5) + 2.25*np.cos(I)**2 * np.cos(0.5*I)**(-4))**(0.5)
	return f_O1(a) * Q_a_inv

#See e.g. Schureman equation 149
def f_Modd(a, n):
	return f_M2(a) ** (n / 2.0)

#Node factors u, see Table 2 of Schureman.

def u_zero(a):
	return 0.0

def u_Mf(a):
	return -2.0 * a['xi'].value

def u_O1(a):
	return 2.0 * a['xi'].value - a['nu'].value

def u_J1(a):
	return -a['nu'].value

def u_OO1(a):
	return -2.0 * a['xi'].value - a['nu'].value

def u_M2(a):
	return 2.0 * a['xi'].value - 2.0 * a['nu'].value

def u_K1(a):
	return -a['nup'].value

#Schureman 214
def u_L2(a):
	I = d2r*a['I'].value
	P = d2r*a['P'].value
	R = r2d*np.arctan(np.sin(2*P)/(1/6.0 * np.tan(0.5*I) **(-2) -np.cos(2*P)))
	return 2.0 * a['xi'].value - 2.0 * a['nu'].value - R

def u_K2(a):
	return -2.0 * a['nupp'].value

#Schureman 202
def u_M1(a):
	I = d2r*a['I'].value
	P = d2r*a['P'].value
	Q = r2d*np.arctan((5*np.cos(I)-1)/(7*np.cos(I)+1)*np.tan(P))
	return a['xi'].value - a['nu'].value + Q

def u_Modd(a, n):
	return n/2.0 * u_M2(a)
