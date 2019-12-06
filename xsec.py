from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import integrate
import lhapdf

# PDF indices: 1 d, 2 u, 21 g 

JPsiMass = 3.096 #Gev
params = [0.3, 0.55, 4.5, 2]    # lambda, kappa, <P_T>, n
"""
def amplitude(x1, x2, P_T, s, params):
    return params[0]**2*params[1]*s*x1*x2/(JPsiMass**2) * np.exp(-params[1]*np.min(P_T**2, params[2])/(JPsiMass**2)) * (1+ np.heaviside(P_T**2 - params[2]**2) * params[1]/params[3] *(P_T**2 - params[2]**2)/(JPsiMass**2))**(-params[3])

def xsec(P_T, s, params, pdf):
    return 1/2s*integrate
"""
ppdf = lhapdf.mkPDF("nCTEQ15_1_1", 0)
#npdf = [lhapdf.mkPDF("nCTEQ15_208_82", ev) for ev in range(33)]


x_vals = np.logspace(-4, 0, 100)
pdf_vals = np.array([ppdf.xfxQ(1, x, 2)-ppdf.xfxQ(-1, x, 2) for x in x_vals])
updf = lambda x: ppdf.xfxQ(2, x, 2)/x-ppdf.xfxQ(-2, x, 2)/x


print(integrate.quad(updf,0.0000,1))

plt.plot(x_vals, pdf_vals)
plt.xscale("log")
plt.show()


