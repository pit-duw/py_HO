from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
from scipy import integrate
sys.path.append('/usr/local/lib/python2.7/site-packages')
import lhapdf

# PDF indices: 1 d, 2 u, 21 g 

JPsiMass = 3.096 #Gev
params = [0.284, 0.55, 4.5, 2]    # lambda, kappa, <P_T>, n

ppdf = lhapdf.mkPDF("nCTEQ15_208_82", 0)
#ppdf = lhapdf.mkPDF("nCTEQ15_208_82/0")

#print((1+ np.heaviside(3.5**2 - params[2]**2) * params[1]/params[3] *(3.5**2 - params[2]**2)/(JPsiMass**2))**(-params[3]))

def amplitude(P_T, params): # Amplitude (without x dependence -> pulled into integrand)
    return params[0]**2*params[1]/(JPsiMass**2) * np.exp(-params[1]*min(P_T**2, params[2]**2)/(JPsiMass**2)) * (1+ np.heaviside(P_T**2-params[2]**2, 0) * params[1]/params[3] *(P_T**2 - params[2]**2)/(JPsiMass**2))**(-params[3])

def integrand(P_T): # Integrand, including potential x dependence of the amplitude
    return lambda x1: ppdf.xfxQ2(21, x1, JPsiMass**2+P_T**2)

def intsquared(P_T, pdf):
    integral = integrate.quad(integrand(P_T), 0,1, limit = 50000)
    return integral[0]**2


PT_vals = np.linspace(0.5, 13.5, 140)
xs_vals = np.array([0.5*389400*amplitude(pt, params)*intsquared(pt, ppdf) for pt in PT_vals])
#print([(PT_vals[i], xs_vals[i]) for i in range(14)])

plt.plot(PT_vals, xs_vals)
plt.yscale("log")
plt.savefig("test.png")
plt.show()