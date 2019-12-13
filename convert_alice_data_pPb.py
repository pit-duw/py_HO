import numpy as np
import yaml

with open("./alice_data/alice_pPb_fwdbwd/data2.yaml") as data_file:
    data = yaml.load(data_file, Loader=yaml.FullLoader)
    pt_mins = [float(data["independent_variables"][0]["values"][i]["low"]) for i in range(len(data["independent_variables"][0]["values"]))]
    pt_maxs = [float(data["independent_variables"][0]["values"][i]["high"]) for i in range(len(data["independent_variables"][0]["values"]))]
    s = float(data["dependent_variables"][0]["qualifiers"][1]["value"])
    y_bounds = [float(y) for y in data["dependent_variables"][0]["qualifiers"][2]["value"].split(" - ")]
    sigma_vals = [float(data["dependent_variables"][0]["values"][i]["value"]) for i in range(len(data["dependent_variables"][0]["values"]))]
    sigma_errs = [np.sqrt(float(data["dependent_variables"][0]["values"][i]["errors"][0]["symerror"])**2+float(data["dependent_variables"][0]["values"][i]["errors"][1]["symerror"])**2+float(data["dependent_variables"][0]["values"][i]["errors"][2]["symerror"])) for i in range(len(data["dependent_variables"][0]["values"]))]
    nbins = len(sigma_vals)
    

def desc(y_bound):
    return "# experimental data d^2sigma/dpT/dy\n# energy ylow yup nbin\n{:.1f} {:.1f} {:.1f} {}\n# ptlow ptup value err\n".format(s, y_bound[0], y_bound[1], nbins)

def line(ptmin, ptmax, sigma, err):
    return "{} {} {} {}\n".format(ptmin, ptmax, sigma, err)


with open("alicepPbdatafwdbwd2.dat", "w+") as file:
    file.write(desc(y_bounds))
    for i in range(nbins):
        file.write(line(pt_mins[i], pt_maxs[i], sigma_vals[i], sigma_errs[i]))

