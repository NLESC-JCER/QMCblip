from ase import Atoms, units
from ase.atoms import Cell

from qmcblip.champ import CHAMP
from qmcblip.flare.quicksim import quicksim, OTFSettings

# Create a thiophene molecule
atoms = Atoms('C4SH4', [(0,0.71495093597,1.31902341514), (0,-0.71495093597,1.31902341514), (0,-1.24055334534,0.05119099544), (0,1.24055334534,0.05119099544), (0,0,-1.15285178278), (0,1.32194923477,2.21441153704), (0,-1.32194923477,2.21441153704), (0,2.27909548764,-0.24288695123), (0,-2.27909548764,-0.24288695123)])
atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
atoms.pbc=[True, True, True]

# Load the FLARE settings and use FLARE++ (and not FLARE)
settings = OTFSettings(theory=OTFSettings.FLAREPP())

# Load the CHAMP calculator
calc = CHAMP(champ_loc="/home/user/bin/vmc.mov1", ncore=4)

# Do a simulation with CHAMP and FLARE++
# With a timestep of 0.5fs and 100 steps
quicksim(atoms, 0.5, 100, calc, settings)
