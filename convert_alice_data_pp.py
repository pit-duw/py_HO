import numpy as np
import yaml

with open("./alice_data/alice_pp/data2.yaml") as data_file:
    data = yaml.load(data_file, Loader=yaml.FullLoader)
    pt_vals = [float(data["independent_variables"][0]["values"][i]["value"]) for i in range(len(data["independent_variables"][0]["values"]))]
    s = float(data["dependent_variables"][0]["qualifiers"][2]["value"])
    y_bounds = [float(y) for y in data["dependent_variables"][0]["qualifiers"][0]["value"].split(" TO ")]
    sigma_vals = [float(data["dependent_variables"][0]["values"][i]["value"]) for i in range(len(data["dependent_variables"][0]["values"]))]
    sigma_errs = [np.sqrt(float(data["dependent_variables"][0]["values"][i]["errors"][0]["symerror"])**2+float(data["dependent_variables"][0]["values"][i]["errors"][1]["symerror"])**2+float(data["dependent_variables"][0]["values"][i]["errors"][2]["symerror"])) for i in range(len(data["dependent_variables"][0]["values"]))]
    nbins = len(sigma_vals)
    

def desc(y_bound):
    return "# experimental data d^2sigma/dpT/dy\n# energy ylow yup nbin\n{:.1f} {:.1f} {:.1f} {}\n# ptlow ptup value err\n".format(s, y_bound[0], y_bound[1], nbins)

def line(ptmin, ptmax, sigma, err):
    return "{} {} {} {}\n".format(ptmin, ptmax, sigma, err)


with open("aliceppdata.dat", "w+") as file:
    file.write(desc(y_bounds))
    pt_bounds = []
    for i in range(nbins):
        if i == 0:
            pt_bounds.append([(3*pt_vals[i]-pt_vals[i+1])/2, (pt_vals[i]+pt_vals[i+1])/2])
        else:
            pt_bounds.append([pt_bounds[i-1][1], 2*pt_vals[i]-pt_bounds[i-1][1]])

    for i in range(nbins):
        file.write(line(pt_bounds[i][0], pt_bounds[i][1], sigma_vals[i], sigma_errs[i]))

