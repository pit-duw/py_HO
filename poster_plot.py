from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import os
from random import shuffle

mpl.rcParams['axes.linewidth'] = 1.2
#plt.style.use("seaborn")
#theory_folder = "./helaconia/hessian/alice_theory/"
#data_folder = "./helaconia/hessian/alice_pp_pPb/"
#theory_folder = "./helaconia/hessian/interpolated_ppfit/"
#data_folder = "./helaconia/hessian/old/"
theory_folder = "./helaconia/hessian/interpolated_alice/"
data_folder = "./helaconia/hessian/y_corr/"

normalization = 1000/208

with open(data_folder+"PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/input/data_list_1.inp", "r") as file:
    data_files = file.read().splitlines()


data_sets = []
data = []
for file in data_files:
    with open(theory_folder+"PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/input/"+file, "r") as f:
        splits = f.read().split("# experimental data d^2sigma/dpT/dy")
        for i in range(1,len(splits)):
            data_sets.append(splits[i].split("\n")[0:3:2])
            data.append(splits[i].split("\n")[4:])


pt_mins = []
pt_maxs = []
for d in data:
    pt_mins.append(np.array([float(line.split(" ")[0]) for line in d[:-1]]))
    pt_maxs.append(np.array([float(line.split(" ")[1]) for line in d[:-1]]))

#pt_mins = np.array(pt_mins)
#pt_maxs = np.array(pt_maxs)

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

colors = [(0.8,0.2,0.5),(0.5,0.8,0.2),(0.2,0.5,0.8)]

s_list = []
y_list = []

figure = plt.figure(1)
for j in range(1,1+len(data_sets)):
    if j == 1:
        continue
    else:
        normalization = 1000/208
    p_array = np.loadtxt(theory_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j+len(data_sets)))[:,0:2]
    pp_data_array = np.loadtxt(data_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,0:2]
    pt_array_p = np.loadtxt(theory_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j+len(data_sets)))[:,0]
    data_array = np.loadtxt(data_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,3:5]
    #pt_array_data = np.loadtxt(data_folder + "PROC_HO_0/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(j))[:,0]
    y_bounds = [data_sets[j-1][1].split(" ")[1], data_sets[j-1][1].split(" ")[2]]
    s = data_sets[j-1][1].split(" ")[0]
    ev_list = []
    s_list.append(s)
    y_list.append(y_bounds)
    
    for i in range(1, 2*N_ev+1)[::-1]:
        array = np.loadtxt(theory_folder + "PROC_HO_{}/P0_addon_fit_pp_psiX_CrystalBall/output/comparison_{}.dat".format(i, j+len(data_sets)))
        ev_list.append(array[:,0:2])
        #plt.plot(array[:,0], array[:,1], ".",  linewidth = 1.5, color= colors[(i-1)//2], zorder = 20)
    
    errors = np.zeros(len(ev_list[0][:,0]))
    for ev in range(0, N_ev*2, 2):
        errors += (ev_list[ev][:,1]-ev_list[ev+1][:,1])**2
        
    errors = np.sqrt(errors)

    frame1 = figure.add_axes((0.15,0.1+1.6/3-(j-2)*0.8/3,0.8,0.8/3))
    
    #plt.title(r"$\sqrt{s} = $" + "{} GeV                 ".format(s) + "{} < $y$ < {}".format(y_bounds[0], y_bounds[1]), fontsize  = 12)
    plt.xlim(0,60)
    plt.plot(pt_array_p, p_array[:,1], "-", zorder = 99, color = (0.0,0.0,0.0))
    plt.errorbar((pt_mins[j-1]+pt_maxs[j-1])/2, data_array[:,0]*normalization, yerr = data_array[:,1]*normalization, xerr = (pt_maxs[j-1]-pt_mins[j-1])/2, zorder = 109, color=colors[j-2], fmt="", capsize=2)
    plt.fill_between(pt_array_p, (ev_list[0][:,1]+ev_list[1][:,1])/2-errors, (ev_list[0][:,1]+ev_list[1][:,1])/2+errors, color = (0.0,0.0,0.0), alpha = 0.3)
    plt.scatter((pt_mins[j-1]+pt_maxs[j-1])/2, data_array[:,0]*normalization, marker="o", zorder = 999, color = colors[j-2], s = 20, label = r"$\sqrt{s} = $" + "{} TeV;   ".format(int(float(s))//1000) + "$y\in[{}, {}]$".format(y_bounds[0], y_bounds[1]))
    plt.yscale("log")
    plt.xlim(0, 20)
    plt.ylim(0.2, 2500)
    plt.legend(frameon=True, loc = "lower left")
    #plt.ylim(min(data_array[:,0])/4, max(data_array[:,0])*4)
    #plt.xticks([])
    plt.ylabel(r"$\frac{d^2\sigma}{dP_Tdy} [\frac{nb}{GeV}]$", fontsize = 12)
    plt.grid(linestyle= "--")
    if j == 4:
        plt.xticks([0,4,8,12,16,20], [0,4,8,12,16,20])
    else:
        plt.xticks([0,4,8,12,16,20], [" ", " ", " ", " ", " ", " "])
    """
    #frame2 = figure.add_axes((0.1,0.1,0.8,0.2))
    #plt.fill_between(pt_array_p, ((ev_list[0][:,1]+ev_list[1][:,1])/2-errors)/((ev_list[0][:,1]+ev_list[1][:,1])/2), ((ev_list[0][:,1]+ev_list[1][:,1])/2+errors)/((ev_list[0][:,1]+ev_list[1][:,1])/2), color = (0.8,0.2,0.5), alpha = 0.3)
    plt.scatter((pt_mins[j-1]+pt_maxs[j-1])/2, data_array[:,0]/pp_data_array[:,1]*normalization, marker="o", zorder = 999, color = (0,0,0), s = 20, label = data_sets[j-1][0][2:-1])
    plt.errorbar((pt_mins[j-1]+pt_maxs[j-1])/2, data_array[:,0]*normalization/pp_data_array[:,1], yerr = data_array[:,1]*normalization/pp_data_array[:,1], xerr = (pt_maxs[j-1]-pt_mins[j-1])/2, zorder = 109, color = (0,0,0), fmt="", capsize=2)
    plt.fill_between(pt_array_p, ((ev_list[0][:,1]+ev_list[1][:,1])/2-errors)/p_array[:,1], ((ev_list[0][:,1]+ev_list[1][:,1])/2+errors)/p_array[:,1], color = (0.8,0.2,0.5), alpha = 0.3)
    plt.ylabel("Data/pp-Theory")
    plt.plot([-1,99999], [1,1], "-", linewidth = "2", color = (0.8,0.2,0.5))
    plt.xlim(0, min(max(pt_maxs[j-1])*1.05,20))
    plt.ylim(0.,2.)
    plt.grid(linestyle= "--")
    """
    #plt.yscale("log")
#plt.scatter((pt_mins[j-1]+pt_maxs[j-1])/2+1000, data_array[:,0]*normalization, marker="o", zorder = 999, color = colors[0], s = 20, label = r"$\sqrt{s} = $" + "{} TeV, ".format(float(s_list[0])//1000) + "{} < $y$ < {}".format(y_list[0][0], y_list[0][1]))
#plt.scatter((pt_mins[j-1]+pt_maxs[j-1])/2+1000, data_array[:,0]*normalization, marker="o", zorder = 999, color = colors[1], s = 20, label = r"$\sqrt{s} = $" + "{} TeV, ".format(float(s_list[1])//1000) + "{} < $y$ < {}".format(y_list[1][1], y_list[1][1]))
plt.xlabel("$P_T$ [GeV]", fontsize = 12)
frame1 = figure.add_axes((0.15,0.9,0.8,0.0001))
plt.plot(pt_array_p, p_array[:,1], "-", zorder = 99, color = (0.0,0.0,0.0), label = "nCTEQ15 proton")
plt.fill_between(pt_array_p, (ev_list[0][:,1]+ev_list[1][:,1])/2-errors, (ev_list[0][:,1]+ev_list[1][:,1])/2+errors, color = (0.0,0.0,0.0), alpha = 0.3, label="Helac-Onia fit uncertainty")
plt.xticks([])
plt.yticks([])
plt.legend()
#plt.legend()
plt.savefig("plots/poster", dpi = 300)
plt.show()