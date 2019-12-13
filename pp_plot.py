from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from random import shuffle

#plt.style.use("seaborn")
#theory_folder = "./helaconia/hessian/alice_theory/"
#data_folder = "./helaconia/hessian/alice_pp_pPb/"
#theory_folder = "./helaconia/hessian/interpolated_ppfit/"
#data_folder = "./helaconia/hessian/old/"
theory_folder = "./helaconia/hessian/alice_pp_pPb/"
data_folder = "./helaconia/hessian/alice_pp_pPb/"

normalization = 1000/208

with open(data_folder+"PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/input/data_list_1.inp", "r") as file:
    data_files = file.read().splitlines()


data_sets = []
for file in data_files:
    with open(data_folder+"PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/input/"+file, "r") as f:
        splits = f.read().split("# experimental data d^2sigma/dpT/dy")
        for i in range(1,len(splits)):
            data_sets.append(splits[i].split("\n")[0:3:2])


#unity = [np.array([0, 100]), np.array([1,1])]
y_bounds = [[2, 2.5],[2.5, 3],[3, 3.5],[3.5, 4],[4, 4.5],[0, 0.75],[0.75, 1.5],[1.5, 2],[2, 2.4],[0, 0.9],[0.9, 1.2],[1.2, 1.6],[1.6, 2.1],[2.1, 2.4]]



lines = {"linestyle": "None"}
plt.rc("lines", **lines)

N_ev = 3
n = int(np.ceil(N_ev**(1/2)))
colors = []
for k in range(n):
    for j in range(n): 
        colors.append( (0.45, j/(n-1)*0.45+0.3, k/(n-1)*0.45+0.3) ) 

shuffle(colors)

for j in range(2,1+len(data_sets)):
    p_array = np.loadtxt(theory_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,0:2]
    pt_array_p = np.loadtxt(theory_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,0]
    data_array = np.loadtxt(data_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,3:5]
    pt_array_data = np.loadtxt(data_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,0]
    y_bounds = [data_sets[j-1][1].split(" ")[1], data_sets[j-1][1].split(" ")[2]]
    s = data_sets[j-1][1].split(" ")[0]
    ev_list = []
    
    for i in range(1, 2*N_ev+1)[::-1]:
        array = np.loadtxt(theory_folder + "PROC_HO_{}/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(i, j))
        ev_list.append(array[:,0:2])
        #plt.plot(array[:,0], array[:,1], ".",  linewidth = 1.5, color= colors[(i-1)//2], zorder = 20)
    
    errors = np.zeros(len(ev_list[0][:,0]))
    for ev in range(0, N_ev*2, 2):
        errors += (ev_list[ev][:,1]-ev_list[ev+1][:,1])**2
        
    errors = np.sqrt(errors)

    figure = plt.figure(1)
    frame1 = figure.add_axes((0.1,0.35,0.8,0.6))
    plt.title(r"$\sqrt{s} = $" + "{} GeV                 ".format(s) + "{} < $y$ < {}".format(y_bounds[0], y_bounds[1]), fontsize  = 12)
    plt.xlim(0,60)
    plt.plot(pt_array_p, p_array[:,1], "-", zorder = 99, color = (0.8,0.2,0.5), label = "nCTEQ15 fit")
    plt.scatter(pt_array_data, data_array[:,0]*normalization, marker="o", zorder = 999, color = (0,0,0), s = 20, label = data_sets[j-1][0][2:-1])
    plt.errorbar(pt_array_data, data_array[:,0]*normalization, yerr = data_array[:,1]*normalization, zorder = 109, color = (0,0,0), fmt="", capsize=2)
    plt.fill_between(pt_array_p, (ev_list[0][:,1]+ev_list[1][:,1])/2-errors, (ev_list[0][:,1]+ev_list[1][:,1])/2+errors, color = (0.8,0.2,0.5), alpha = 0.3, label="Uncertainty from HO fit")
    plt.yscale("log")
    plt.legend()
    plt.xlim(0, max(pt_array_data)*1.05)
    #plt.ylim(min(data_array[:,0])/4, max(data_array[:,0])*4)
    #plt.xticks([])
    plt.ylabel(r"$\frac{d^2\sigma}{dP_Tdy} [\frac{nb}{GeV}]$", fontsize = 12)
    plt.grid(linestyle= "--")

    frame2 = figure.add_axes((0.1,0.1,0.8,0.2))
    #plt.fill_between(pt_array_p, ((ev_list[0][:,1]+ev_list[1][:,1])/2-errors)/((ev_list[0][:,1]+ev_list[1][:,1])/2), ((ev_list[0][:,1]+ev_list[1][:,1])/2+errors)/((ev_list[0][:,1]+ev_list[1][:,1])/2), color = (0.8,0.2,0.5), alpha = 0.3)
    plt.scatter(pt_array_data, data_array[:,0]/p_array[:,1]*normalization, marker="o", zorder = 999, color = (0,0,0), s = 20, label = data_sets[j-1][0][2:-1])
    plt.errorbar(pt_array_data, data_array[:,0]*normalization/p_array[:,1], yerr = data_array[:,1]*normalization/p_array[:,1], zorder = 109, color = (0,0,0), fmt="", capsize=2)
    plt.fill_between(pt_array_p, ((ev_list[0][:,1]+ev_list[1][:,1])/2-errors)/p_array[:,1], ((ev_list[0][:,1]+ev_list[1][:,1])/2+errors)/p_array[:,1], color = (0.8,0.2,0.5), alpha = 0.3)
    plt.ylabel("Data/pp-Theory")
    plt.plot([-1,99999], [1,1], "-", linewidth = "2", color = (0,0,0))
    plt.xlabel("$P_T$ [GeV]", fontsize = 12)
    plt.xlim(0, max(pt_array_data)*1.05)
    plt.ylim(0.,2.)
    plt.grid(linestyle= "--")
    #plt.yscale("log")
    #plt.savefig("plots/alice_pPb_comparison_{}".format(j), dpi = 300)
    plt.show()