CHAMP Calculator
-----------------

Setting up the calculator
^^^^^^^^^^^^^^^^^^^^^^^^^

The CHAMP calculator can used in ASE for any kind of simulation. 
After importing the calculator, a new CHAMP calculator object can be made using:

>>> from qmcblip.champ import CHAMP
>>> champ_calc = CHAMP(champ_loc="/home/user/champ/bin/vmc.mov1", force_file="write_forces", ncore=4)

To view all the parameters that can be passed to the calculator, see the :class:`CHAMP <qmcblip.champ.CHAMP>` class documentation.

To attach this calculator to a molecule and do Molecular Dynamics, we can do:

>>> from ase.md.verlet import VelocityVerlet
>>> from ase.build import molecule
>>>
>>> # Make a H2 molecule
>>> atoms = molecule("H2")
>>> 
>>> # Attach the calculator
>>> atoms.calc = champ_calc
>>> 
>>> # Setup and run MD
>>> dyn = VelocityVerlet(atoms, 0.5 * units.fs) 
>>> dyn.run(30)


Configuring CHAMP
^^^^^^^^^^^^^^^^^

There are two methods of configuring CHAMP. You can have a ``vmc.inp`` file in the same directory as the Python script.
This ``vmc.inp`` file should contain all the neccessery settings. 
The second method is by passing a :obj:`Settings<qmcblip.champio.Settings>` object to the CHAMP calculator, using the ``settings`` keyword argument.
A :obj:`Settings<qmcblip.champio.Settings>` object can be created by reading from a ``vmc.inp`` file (see :func:`qmcblip.champio.Settings.read`):

>>> from qmcblip.champio import Settings
>>> settings = Settings.read('vmc.inp')
>>>
>>> champ_calc = CHAMP(settings=settings)

It can also be created by creating it explicitly in Python and passing the neccessery arguments:

>>> settings = Settings(...)

See the :class:`Settings<qmcblip.champio.Settings>` documentation for an overview of which parameters are required.
Once you have a :obj:`Settings<qmcblip.champio.Settings>` object, it can also be printed to a ``vmc.inp`` file (see :func:`qmcblip.champio.Settings.write`):

>>> settings.write('vmc.inp')
