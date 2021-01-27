#! /usr/bin/env python
import numpy as np
import subprocess
import os
from numpy import linalg as LA

central_value_pdf = 100000

def get_central_values(file):
    with open(file, "r") as file:
        settings = file.read()

    settings_lines = settings.splitlines()
    params = []
    for i in range(4):
        params.append(float(settings_lines[i].split(" ")[0][0:-2]))

    return params


def generate_fit_param_card(params, file):
    with open(file, "w") as file:
        file.write(str(params[0])+"d0  0.0d0 0d0 10d0 # kappa\n")
        file.write(str(params[1])+"d0  0.0d0 0d0 10d0 # lam\n")
        file.write(str(params[2])+"d0  0.0d0 0d0 20d0 # <pt>\n")
        file.write(str(params[3])+"d0  0.0d0 0d0 10d0 # n\n")
        file.write("# initial-value step-size lower-bound upper-bound\n# step-size=0d0 means they are fixed and the bounds irrelevant\n# lower and upper bound of MINUIT parameters\n# if both of them are 0d0, there is no limitation")

def calc_results():
    p = subprocess.Popen("../helaconia/cluster/bin/ho_cluster", 
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)

    p.stdin.write("generate addon 3\n")
    p.stdin.write("launch")
    p.communicate()

def get_chisq(params):
    print("Calculating chisq at parameters: ", params)
    generate_fit_param_card(params, "../helaconia/addon/fit_pp_psiX_CrystalBall/input/fit_param_card.inp")
    calc_results()
    folders = os.listdir(".")
    latest = 0
    for folder in folders:
        if "PROC_HO_" in folder:
            latest = max(latest, int(folder.replace("PROC_HO_", "")))
    chisq = 0.0
    with open ("PROC_HO_"+str(latest)+"/P0_addon_fit_pp_psiX_CrystalBall/output/fit_results.out") as out:
        output = out.read()
        chisq= float(output.splitlines()[-1].split()[-1])

    print("chisq: ", chisq)
    generate_fit_param_card(central, "../helaconia/addon/fit_pp_psiX_CrystalBall/input/fit_param_card.inp")
    return chisq

def unit_vec(i):
    uv = np.zeros(4)
    uv[i] = 1.0
    return uv

def get_hessian(central):
    hess = np.zeros((3,3))
    dp = 0.0005
    cv = np.array(central)
    for i in range(3):
        hess[i,i] = (get_chisq(cv+unit_vec(i)*dp) - 2*chisq_c + get_chisq(cv-unit_vec(i)*dp))/(dp**2)
        for j in range(i):
            if not i == j:
                hess[i,j] = (get_chisq(cv+unit_vec(i)*dp+unit_vec(j)*dp) + get_chisq(cv-unit_vec(i)*dp-unit_vec(j)*dp) - get_chisq(cv+unit_vec(j)*dp-unit_vec(i)*dp) - get_chisq(cv+unit_vec(i)*dp-unit_vec(j)*dp))/(4*dp**2)
                hess[j,i] = hess[i,j]
    np.savetxt("hessian.txt", hess)
    return hess


central = get_central_values("./backup_fit_param_card.inp")
print("Central values of parameters: ", central)
#generate_fit_param_card(central, "./backup_fit_param_card.inp")
get_chisq(central)
chisq_c = 177.24
T = np.sqrt(np.sqrt(2*chisq_c))
t = T/np.sqrt(3)
print("chisq at central values: ", chisq_c)
#print(get_hessian(central))
    
eigenvals, eigenvecs = LA.eig(np.loadtxt("hessian.txt"))

#print(eigenvals, eigenvecs)

n_eigenvecs = np.zeros((3,4))
for i in range(3):
    n_eigenvecs[i, 0:3] = eigenvecs[i]*np.sqrt(np.abs(2/eigenvals[i]))
    if i == 0:
        n_eigenvecs[i] = n_eigenvecs[i] * 1.1
    if i == 1:
        n_eigenvecs[i] = n_eigenvecs[i] * 1.4
    if i == 2:
        n_eigenvecs[i] = n_eigenvecs[i] /35

print(n_eigenvecs[0])
print(n_eigenvecs[1])
print(n_eigenvecs[2])

#np.savetxt("eigenvecs.txt", n_eigenvecs)



for i in range(3):
    print("dchisq of eigenvector " + str(i) + " in + direction: " + str(get_chisq(central+t*n_eigenvecs[i])))
    print("dchisq of eigenvector " + str(i) + " in - direction: " + str(get_chisq(central-t*n_eigenvecs[i])))

