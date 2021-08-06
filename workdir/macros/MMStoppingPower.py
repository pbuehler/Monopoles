import numpy as np
from scipy.integrate import quad

def Im (Z):       # magnetic mean ionization potential (eV)
  val = 0;
  if Z <= 13:
    val = 12.*Z + 7.
  else:
    val = 0.76*Z + 58.8*np.power(Z,-0.19)
    
  return val

def K (g):
  val = 0.406
  if g > 1:
    val = 0.346
  
  return val
  
def B (g):
  val = 0.248
  if g > 1:
    val = 0.627
  if g > 2:
    val = 1.002
  if g > 3:
    val = 1.243
  if g > 4:
    val =1.464
    
  return val
  
def dm (beta, Z):   # density-effect correction
  val = 0.
  
  return val

# g: gD [68.5 e]
# beta: v/c [1]
#
def dEdx (g, beta):
  
  # some material parameters
  me =    # mass of electron [GeV/c^2]
  N =     # electron density
  e = 0   # electric charge
  
  gamma = 1./np.sqrt(1.-beta*beta)
  bg = beta*gamma
  
  v0 = 4*np.pi*N*g*g*e*e/me
  v1 = 0.5*np.log(2.*me*bg*bg/(2.*np.power(Im(Z),2.))
  v1 += K(g)/2.-0.5-dm(beta, Z)/2.-B(g)
  
  return v0*v1

a = 2
b = 1
I = quad(integrand, 0, 1, args=(a,b))
I
(1.6666666666666667, 1.8503717077085944e-14)
