from os.path import exists
from os import remove, rename
from glob import glob

wf = ["determinants", "orbitals", "jastrow", "jastrow_der", "molecule", "basis_num_info", "symmetry"]

def check_tags(tags):
    # Check if the argument is a dict
    if (type(tags) is not dict):
        raise TypeError("The input is not a dictionary!")

    # Check if the general module is present
    if ("general" not in tags) or (type(tags['general']) is not dict):
        raise ValueError("Make sure your input contains the 'general' module!")

    # Check if the wavefunction and basis files are present
    for item in wf:
        if (item not in tags) and (item != "symmetry"):
            raise ValueError("Input does not contain " + item + "!")
        if (not exists(tags[item])) and (item != "molecule"):
            raise FileNotFoundError(tags[item] + " was not found!")


def write_input(tags, filename="vmc.inp"):
    # Check the validity of the input
    check_tags(tags)

    f = open(filename, 'w')

    for key, value in tags.items():
        if (type(value) is dict):
            f.write("\n%module " + key + "\n")
            for key2, value2 in value.items():
                f.write("\t" + key2 + " " + str(value2) + "\n")
            f.write("%endmodule\n\n")
        elif (key in wf):
            f.write("load " + key + " " + str(value) + "\n")
        else:
            f.write(key + " " + str(value) + "\n")


    f.close()

def read_input(filename="vmc.inp"):
    # Check if the file exists
    if not exists(filename):
        raise FileNotFoundError(filename + " was not found!")

    output = dict()

    f = open(filename, 'r')

    curmod = ""
    inmod = False

    for line in f:
        if not line.strip():
            continue
        elif line.startswith("%module"):
            curmod = line.removeprefix("%module ").removesuffix('\n')
            inmod = True
            output[curmod] = dict()
        elif line.startswith("%endmodule"):
            inmod = False
        elif line.startswith("load"):
            temp = line.removeprefix("load ").removesuffix('\n')
            temp = temp.split()
            key = temp[0]
            temp.pop(0)
            value = ' '.join(temp)
            output[key] = value
        else:
            temp = line.removesuffix('\n').split()
            key = temp[0]
            temp.pop(0)
            value = ' '.join(temp)
            if inmod:
                output[curmod][key] = value
            else:
                output[key] = value

    f.close()

    return output

def cleanup(*args):
    standard = ["force_analytic", "restart_vmc", "output.log", "parser.log", "orbitals_optimal.1.*",
                "jastrow_optimal.1.*", "det_optimal.1.*", "geo_optimal.*", "geo_optimal_final",
                "mc_configs_start", "nohup.out", "vmc_temp.inp"]
    for item in args:
        if item in standard:
            standard.remove(item)
        else:
            standard.append(item)

    for file in standard:
        if (file.find('*') != -1):
            for f in glob(file):
                remove(f)
        else:
            if exists(file):
                remove(file)

def use_opt_wf(filename="vmc.inp"):
    opt = ["det_optimal.1.iter*", "orbitals_optimal.1.iter*", "jastrow_optimal.1.iter*"]
    keys = ["determinants", "orbitals", "jastrow"]

    tags = read_input(filename)

    for ind, name in enumerate(opt):
        num = len(glob(name))
        if num > 0:
            opt[ind] = opt[ind].strip('*') + str(num)
            tags[keys[ind]] = opt[ind]

    write_input(tags, filename)