from __future__ import division
import numpy as np
#import matplotlib.pyplot as plt
import os
import sys
from scipy import integrate
sys.path.append('/usr/local/lib/python2.7/site-packages')
import lhapdf

# PDF indices: 1 d, 2 u, 21 g 

JPsiMass = 3.096 #Gev
params = [0.55, 0.283, 4.5, 2]    # lambda, kappa, <P_T>, n

ppdf = lhapdf.mkPDF("nCTEQ15_208_82", 0)

#print((1+ np.heaviside(3.5**2 - params[2]**2) * params[1]/params[3] *(3.5**2 - params[2]**2)/(JPsiMass**2))**(-params[3]))

def amplitude(x1, x2, P_T, s, params):
    return params[0]**2*params[1]*s*x1*x2/(JPsiMass**2) * np.exp(-params[1]*min(P_T**2, params[2])/(JPsiMass**2)) * (1+ np.heaviside(P_T**2, params[2]**2) * params[1]/params[3] *(P_T**2 - params[2]**2)/(JPsiMass**2))**(-params[3])

def integrand(P_T, s, params):
    return lambda x1, x2: ppdf.xfxQ2(21, x1, JPsiMass**2+P_T**2)*ppdf.xfxQ2(21, x2, JPsiMass**2+P_T**2)*amplitude(x1, x2, P_T, s, params)

def xsec(P_T, s, params, pdf):
    return 1/(2*s)*integrate.nquad(integrand(P_T, s, params), [[0,1],[0,1]])[0]


#npdf = [lhapdf.mkPDF("nCTEQ15_208_82", ev) for ev in range(33)]
print(389.4*xsec(3, 7000**2, params, ppdf))

1/0

x_vals = np.logspace(-4, 0, 100)
pdf_vals = np.array([ppdf.xfxQ(1, x, 2)-ppdf.xfxQ(-1, x, 2) for x in x_vals])
updf = lambda x: ppdf.xfxQ(2, x, 2)/x-ppdf.xfxQ(-2, x, 2)/x


print(integrate.quad(updf,0.0000,1))

#plt.plot(x_vals, pdf_vals)
#plt.xscale("log")
#plt.show()


