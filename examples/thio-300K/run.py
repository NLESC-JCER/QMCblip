# Import numpy and matplotlib
import numpy as np
from numpy.random import random

# flare++ imports
import flare_pp._C_flare as flare_pp
from flare_pp.sparse_gp import SGP_Wrapper
from flare_pp.sparse_gp_calculator import SGP_Calculator

# flare imports
from qmcblip.otf import C_ASE_OTF as ASE_OTF
from flare import otf_parser

# ASE imports
from ase import Atoms, units
from ase.atoms import Cell
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, \
    Stationary, ZeroRotation

from qmcblip.champ import CHAMP

np.random.seed(12345)

#atoms = molecule("H2")
atoms = Atoms('C4SH4', [(0,0.71495093597,1.31902341514), (0,-0.71495093597,1.31902341514), (0,-1.24055334534,0.05119099544), (0,1.24055334534,0.05119099544), (0,0,-1.15285178278), (0,1.32194923477,2.21441153704), (0,-1.32194923477,2.21441153704), (0,2.27909548764,-0.24288695123), (0,-2.27909548764,-0.24288695123)])
atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
atoms.pbc=[True, True, True]
n_atoms = len(atoms)

# Set MD parameters.
md_engine = "CustomVerlet"
md_dict = {}

MaxwellBoltzmannDistribution(atoms, temperature_K=300)
Stationary(atoms)  # zero linear momentum
ZeroRotation(atoms)  # zero angular momentum

# Create sparse GP model.
species_map = {16: 0, 6: 1, 1: 2} 
cutoff = 5.0  # in A
sigma = 2.0  # in eV
power = 2  # power of the dot product kernel
kernel = flare_pp.NormalizedDotProduct(sigma, power)
cutoff_function = "quadratic"
many_body_cutoffs = [cutoff]
radial_basis = "chebyshev"
radial_hyps = [0., cutoff]
cutoff_hyps = []
n_species = 3
N = 12
lmax = 3
descriptor_settings = [n_species, N, lmax]
descriptor_calculator = flare_pp.B2(
  radial_basis,
  cutoff_function,
  radial_hyps,
  cutoff_hyps,
  descriptor_settings
)

# Set the noise values.
sigma_e = 0.12 * n_atoms  # Energy noise (in kcal/mol, so about 5 meV/atom)
sigma_f = 0.115  # Force noise (in kcal/mol/A, so about 5 meV/A)
sigma_s = 0.014  # Stress noise (in kcal/A^3, so about 0.1 GPa)
# Choose uncertainty type.
# Other options are "DTC" (Deterministic Training Conditional) or
# "SOR" (Subset of Regressors).
variance_type = "local"  # Compute uncertainties on local energies (normalized)

# Choose settings for hyperparameter optimization.
max_iterations = 10  # Max number of BFGS iterations during optimization
opt_method = "L-BFGS-B"  # Method used for hyperparameter optimization

# Bounds for hyperparameter optimization.
# Keeps the energy noise from going to zero.
bounds = [(None, None), (sigma_e, None), (None, None), (None, None)]

# Create a model wrapper that is compatible with the flare code.
gp_model = SGP_Wrapper([kernel], [descriptor_calculator], cutoff,
                    sigma_e, sigma_f, sigma_s, species_map,
                    variance_type=variance_type,
                    stress_training=False,
                    opt_method=opt_method,
                    bounds=bounds,
                    max_iterations=max_iterations)

# Create an ASE calculator based on the GP model.
flare_calculator = SGP_Calculator(gp_model)

# Set up OTF object.
init_atoms = list(range(n_atoms))  # Initial environments to include in the sparse set
output_name = 'thio'  # Name of the output file
std_tolerance_factor = -0.02  # Uncertainty tolerance for calling DFT
freeze_hyps = 10  # Freeze hyperparameter optimization after this many DFT calls
min_steps_with_model = 0  # Min number of steps between DFT calls
update_style = "threshold"  # Strategy for adding sparse environments
update_threshold = 0.005  # Threshold for determining which sparse environments to add

otf_params = {'init_atoms': init_atoms, 'output_name': output_name,
              'std_tolerance_factor': std_tolerance_factor,
              'freeze_hyps': freeze_hyps,
              'force_only': False,
              'min_steps_with_model': min_steps_with_model,
              'update_style': update_style,
              'update_threshold': update_threshold}
              
champ_calc = CHAMP(champ_loc="/home/user/bin/vmc.mov1", force_file="write_forces", ncore=60, nodefile="nodefile")

changes = []

# Create OTF object.
timestep = 0.5 * units.fs
number_of_steps = 1000
test_otf = ASE_OTF(
    atoms, timestep=timestep, number_of_steps=number_of_steps,
    dft_calc=champ_calc, md_engine=md_engine, md_kwargs=md_dict,
    update_settings=changes, calculator=flare_calculator, **otf_params)

# Run on-the-fly dynamics.
test_otf.run()
