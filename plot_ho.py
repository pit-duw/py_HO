from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from random import shuffle

#plt.style.use("seaborn")

unity = [np.array([0, 100]), np.array([1,1])]
y_bounds = [[2, 2.5],[2.5, 3],[3, 3.5],[3.5, 4],[4, 4.5],[0, 0.75],[0.75, 1.5],[1.5, 2],[2, 2.4],[0, 0.9],[0.9, 1.2],[1.2, 1.6],[1.6, 2.1],[2.1, 2.4]]
#data  = np.loadtxt("./PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_1.dat".format(i))

lines = {"linestyle": "None"}
plt.rc("lines", **lines)

N_ev = 16
n = int(np.ceil(N_ev**(1/3)))
colors = []
for k in range(n+1):
    for j in range(n+1):
        for i in range(n+1):
            colors.append( (i/n*0.45+0.3, j/n*0.45+0.3, k/n*0.45+0.3) ) 

shuffle(colors)

for j in range(1,15):
    p_array = np.loadtxt("./PROC_HO_P/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,1]
    ev_list = []

    for i in range(1, 2*N_ev+1)[::-1]:
        array = np.loadtxt("./PROC_HO_{}/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(i, j))
        
        ev_list.append(array[:,0:2])

        lw = 5 if not i else 1.5
        zo = 100 if not i else 10
        plt.plot(array[:,0], array[:,1]/p_array, ".",  linewidth = lw, color= colors[(i-1)//2], zorder = zo)
        #plt.yscale("log")

    errors = np.zeros(len(ev_list[0][:,0]))
    for ev in range(0, 32, 2):
        errors += (ev_list[ev][:,1]-ev_list[ev+1][:,1])**2
        
    errors = np.sqrt(errors)

    plt.xlabel(r"$P_T$ [GeV]", fontsize = 12)
    plt.ylabel(r"$R_{pPb} = \left(\frac{d^2\sigma_{pPb}}{dP_Tdy}\right)\cdot \left(\frac{d^2\sigma_{pp}}{dP_Tdy}\right)^{^{-1}}$", fontsize  = 12)
    plt.title(r"$\sqrt{s} = 7$ TeV                 " + "{} < $y$ < {}".format(y_bounds[j-1][0], y_bounds[j-1][1]), fontsize  = 12)
    plt.xlim(0,60)
    plt.plot(unity[0], unity[1], "--", color = (0.5,0.5,0.5), zorder = 1)
    plt.errorbar(ev_list[0][:,0], ev_list[0][:,1]/p_array, yerr = errors/p_array, zorder = 9, color = (0,0,0), fmt="", capsize=2)
    plt.scatter(ev_list[0][:,0], ev_list[0][:,1]/p_array, marker="o", zorder = 999, color = (0,0,0), s = 40)
    plt.savefig("plots/comparison_{}".format(j))
    plt.show()