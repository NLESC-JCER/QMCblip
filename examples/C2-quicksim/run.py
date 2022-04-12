from ase import Atoms, units
from ase.atoms import Cell

from qmcblip.champ import CHAMP
from qmcblip.flare.quicksim import quicksim, OTFSettings

# Create a C2 molecule
atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
atoms.pbc=[True, True, True]

# Load the FLARE settings and use FLARE (and not FLARE++)
settings = OTFSettings(theory=OTFSettings.FLARE())

# Load the CHAMP calculator
calc = CHAMP(champ_loc="/home/user/bin/vmc.mov1", ncore=4)

# Do a simulation with CHAMP and FLARE
# With a timestep of 0.5ps and 100 steps
quicksim(atoms, 0.5, 100, calc, settings)
