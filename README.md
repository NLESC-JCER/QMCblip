# CHAMP-ASE

This repository contains a [ASE](https://gitlab.com/ase/ase) calculator for CHAMP. In this README we will go over the process of using this calculator.

# Requirements

Simply install the following packages using pip (see [this](https://gitlab.com/ase/ase) for the ASE requirements):
```
pip3 install ase numpy
```

## Installing CHAMP

To use this software, you need a special version of CHAMP which can export the forces and energy to a file. For this you need to [ase-coupling](https://github.com/filippi-claudia/champ/tree/ase-coupling) branch. Simply clone this and build in the usual way.

## ASE CHAMP calculator

### Installation
To install this calculator, simply place it in the same folder as your python script (we will later push it to ASE). To import:
```
from champ import CHAMP
```
To use this calculator if the python file is located in another folder (recommended), simply use:
```
import sys
sys.path.append('/path/to/folder/')
from champ import CHAMP
```
### Running
After importing the calculator, a new CHAMP calculator object can be made using:
```
champ_calc = CHAMP(champ_loc="/home/user/champ/bin/vmc.mov1", force_file="write_forces", ncore=4)
```
The CHAMP specific parameters that can be set are:
- `vmc_in` -- The input file for CHAMP (default 'vmc.inp')
- `vmc_out` -- The output file of CHAMP (default 'vmc.out')
- `force_file` -- The file that CHAMP writes the forces and energies to (default 'champ_forces')[^1]
- `pos_file` -- The file from which CHAMP will read to location of the atoms (default 'molecule.xyz')[^2]
- `champ_loc` -- Location of the CHAMP executable (default '/usr/bin/vmc.mov1')
- `nodefile` -- If set, the calculator will run on multiple nodes for CHAMP
- `ncore` -- Amount of cores to run CHAMP on (default 1)

[^1]: As of now there is no way to set this in CHAMP, so 'write_forces' should be used.
[^2]: Should match the name in the CHAMP input file.

To attach this calculator to a molecule and do molecular dynamics, we can do:
```
from ase.md.verlet import VelocityVerlet
from ase.build import molecule

# Make a H2 molecule
atoms = molecule("H2")

# Attach the calculator
atoms.calc = champ_calc

# Setup and run MD
dyn = VelocityVerlet(atoms, 0.5 * units.fs) 
dyn.run(30)
```
