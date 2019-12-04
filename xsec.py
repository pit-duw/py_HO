from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import integrate
import lhapdf

JPsiMass = 3.096 #Gev


ppdf = lhapdf.mkPDF("nCTEQ15_1_1", 0)
npdf = [lhapdf.mkPDF("nCTEQ15_208_82", ev) for ev in range(33)]


x_vals = np.logspace(-4, 0, 100)
pdf_vals = np.array([pdf.xfxQ(1, x, 2)-pdf.xfxQ(-1, x, 2) for x in x_vals])

plt.plot(x_vals, pdf_vals)
plt.xscale("log")
plt.show()


