#! /usr/bin/env python

import subprocess
import os

central_value_pdf = 101100
min_ev_pdf = 101101
max_ev_pdf = 101132

try:
    delete  = input("Write 1 to delete all old PROC folders, anything else to continue without deleting.\n")
except SyntaxError:
    delete = None

if delete == 1:
    print("Deleting old results")
    del_result = os.system('rm -r PROC*')
else: 
    print("Keeping old results.")

def set_pdf(number):
    with open("./helaconia/input/user.inp", "r") as file:
        settings = file.read()

    settings_lines = settings.splitlines()

    for i, line in enumerate(settings_lines):
        if line.startswith("pdf "):
            settings_lines[i] = "pdf " + str(number)

    with open("./helaconia/input/user.inp", "w") as file:
        file.write("\n".join(settings_lines))

def calc_results():
    p = subprocess.Popen("./helaconia/cluster/bin/ho_cluster", 
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)

    p.stdin.write("generate addon 3\n")
    p.stdin.write("launch")
    p.communicate()

set_pdf(central_value_pdf)
calc_results()
print("Calculation for central value done.")

for i in range(min_ev_pdf, max_ev_pdf+1):
    set_pdf(i)
    calc_results()
    print("Calculation for eigenvector" + str(i) +" done.")

set_pdf(central_value_pdf)

os.system('rm ./helaconia/addon/fit_pp_psiX_CrystalBall/input/fit_param_card.inp')
