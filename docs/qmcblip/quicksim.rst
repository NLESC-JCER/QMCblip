Quickly setup
-------------

QMCblip comes with a quicksim function to quickly perform QMC ML FF simulations. To do a quick FLARE simulation, import:

>>> from qmcblip.flare.quicksim import quicksim, OTFSettings

Depending on your system, you might want to use FLARE or FLARE++. To switch between the two, make an :obj:`OTFSettings<qmcblip.flare.quicksim.OTFSettings>` object:

>>> settings = OTFSettings(theory=OTFSettings.FLARE())

In the case you want to do a FLARE++ simulation, do:

>>> settings = OTFSettings(theory=OTFSettings.FLAREPP())

You also need to provide the CHAMP calculator

>>> from qmcblip.champ import CHAMP
>>> calc = CHAMP(champ_loc="/home/user/bin/vmc.mov1", ncore=4)

Finally, you also need to setup the molecular system (remember to put a box around it for FLARE):

>>> from ase import Atoms, units
>>> from ase.atoms import Cell
>>>
>>> atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
>>> atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
>>> atoms.pbc=[True, True, True]

Finally you are ready to perform the simulation. In this case we are performing a simulation of 100 steps with 0.5fs timestep.

>>> quicksim(atoms, 0.5, 100, calc, settings)

You can also supply an extra argument to change the CHAMP settings during the simulation. To change the amount of optimization steps at the 10th MD timestep:

>>> changes = [(10, {'optwf': {'nopt_iter': 10}})]
>>> quicksim(atoms, 0.5, 100, calc, settings, changes=changes)

For more examples, see :doc:`../examples/examples`.


