from os.path import exists
from os import remove, rename
from glob import glob

wf = ["determinants", "orbitals", "jastrow", "jastrow_der", "molecule", "basis_num_info", "symmetry"]

def check_tags(tags):
    """
    Check if a dictionary will satisfy basic neccessities for CHAMP input file.
    This does not check everything.

    Arguments:
    tags -- dictionary containing CHAMP input
    """
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
    """
    Write a dictionary containing the CHAMP configuration to an input file.

    Arguments:
    tags -- dictionary containing CHAMP input. Should satisfy check_tags() function
    filename -- input file to write to
    """
    # Check the validity of the input
    check_tags(tags)

    f = open(filename, 'w')

    for key, value in tags.items():
        # Writing a module
        if (type(value) is dict):
            f.write("\n%module " + key + "\n")
            for key2, value2 in value.items():
                f.write("\t" + key2 + " " + str(value2) + "\n")
            f.write("%endmodule\n\n")
        # Loading in a file
        elif (key in wf):
            f.write("load " + key + " " + str(value) + "\n")
        else:
            f.write(key + " " + str(value) + "\n")

    f.close()

def read_input(filename="vmc.inp"):
    """
    Read the CHAMP input file and convert it to a dictionary format

    Arguments:
    filename -- file of the input file to read and write to

    Output:
    Dictionary containing the CHAMP configuration
    """
    # Check if the file exists
    if not exists(filename):
        raise FileNotFoundError(filename + " was not found!")

    output = dict()

    f = open(filename, 'r')

    # Set some variables to keep track of modules
    curmod = ""
    inmod = False

    for line in f:
        # We skip empty lines
        if not line.strip():
            continue
        # Beginning of a module
        elif line.startswith("%module"):
            curmod = line.removeprefix("%module ").removesuffix('\n')
            inmod = True
            output[curmod] = dict()
        # End of a module
        elif line.startswith("%endmodule"):
            inmod = False
        # Loading in a file
        elif line.startswith("load"):
            temp = line.removeprefix("load ").removesuffix('\n')
            temp = temp.split()
            key = temp[0]
            temp.pop(0)
            value = ' '.join(temp)
            output[key] = value
        # All other tags
        else:
            temp = line.removesuffix('\n').split()
            key = temp[0]
            temp.pop(0)
            value = ' '.join(temp)
            if inmod:
                # If we are in a module, add the tag to that subdictionary.
                output[curmod][key] = value
            else:
                output[key] = value

    f.close()

    # Return the dictionary
    return output

def cleanup(*args):
    """
    Remove files created by CHAMP. Add files as arguments to include (if not in default list)
    or exclude (if in default list) files.
    """
    # Default list of the files that will be removed
    standard = ["force_analytic", "restart_vmc", "output.log", "parser.log", "orbitals_optimal.1.*",
                "jastrow_optimal.1.*", "det_optimal.1.*", "geo_optimal.*", "geo_optimal_final",
                "mc_configs_start", "nohup.out", "vmc_temp.inp"]
    
    # Add or remove file depending on the input arguments
    for item in args:
        if item in standard:
            standard.remove(item)
        else:
            standard.append(item)

    for file in standard:
        # Take care of the wildcard files
        if (file.find('*') != -1):
            for f in glob(file):
                remove(f)
        else:
            if exists(file):
                remove(file)

def use_opt_wf(filename="vmc.inp"):
    """
    Function to replace the orbitals, determinants and jastrow with the optimized files.

    Arguments:
    filename -- file of the input file to read and write to
    """
    opt = ["det_optimal.1.iter*", "orbitals_optimal.1.iter*", "jastrow_optimal.1.iter*"]
    keys = ["determinants", "orbitals", "jastrow"]

    tags = read_input(filename)

    for ind, name in enumerate(opt):
        num = len(glob(name))
        # If there are no optimized WF files, we do not use them
        if num > 0:
            opt[ind] = opt[ind].strip('*') + str(num)
            tags[keys[ind]] = opt[ind]

    write_input(tags, filename)