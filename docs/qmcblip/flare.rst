.. _FLARE: https://github.com/mir-group/flare
.. _`FLARE++`: https://github.com/mir-group/flare_pp

.. _`FLARE documentation`: https://flare.readthedocs.io/en/latest/
.. _`FLARE++ documentation`: https://mir-group.github.io/flare_pp/

FLARE
-----

Introduction
^^^^^^^^^^^^

FLARE_ is an on-the-fly Force Field Machine Learning packages. 
It uses Gaussian Process regression to train the Force Field. 
However, unfortunatly for use it only works for atoms in a box, since it also determines the stress. 
Without a box, you cannot calculate stress. 
We circumvent this by placing the molecule in a big box and letting the CHAMP calculator pass a null tensor. 
second problem is that it only works up to 3-body interactions. 
However for MD we idealy want atleast 4-body interactions. 
This is where `FLARE++`_ comes in. 
FLARE++ is an extension to FLARE to allow for n-body interactions. 
This package contians a modified version of the on-the-fly learning code from FLARE, to make it work well with QMC. 
To this aim, there is also a modified Velocity Verlet scheme which must be used.

Using FLARE
^^^^^^^^^^^
See the `FLARE documentation`_ and `FLARE++ documentation`_ on how to use FLARE together with ASE.
Here we will quickly explain the specific steps that are needed in order to function with CHAMP.

Instead of import the :obj:`ASE_OTF <flare.ase.otf.ASE_OTF>` from FLARE we need to import the :obj:`C_ASE_OTF <qmcblip.otf.C_ASE_OTF>` object:

>>> from qmcblip.otf import C_ASE_OTF as ASE_OTF

Because FLARE expect our molecule to be in a box, we also have to put a box around our molecule:

>>> from ase.build import molecule
>>> from ase.atoms import Cell
>>>
>>> # Make a H2 molecule
>>> atoms = molecule("H2")
>>> atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
>>> atoms.pbc=[True, True, True]

We also use a custom MD engine for our QMC-MD simulations:

>>> md_engine = 'CustomVerlet'
>>> md_kwargs = {}

Additionally, we pass an array which tells the OTF runner if and when we want to update CHAMP settings. 
This array can be empty.

>>> changes = [(10, {'optwf': {'nopt_iter': 2}})]

>>> test_otf = ASE_OTF(atoms, 
>>>                    timestep = 0.5 * units.fs,
>>>                    number_of_steps = 500,
>>>                    dft_calc = champ_calc,
>>>                    md_engine = md_engine,
>>>                    md_kwargs = md_kwargs,
>>>                    update_settings = changes,
>>>                    calculator=flare_calculator,
>>>                    **otf_params)

See the `FLARE documentation`_ for more information about the other parts in the code on this page.

Analyzing results
^^^^^^^^^^^^^^^^^

Due to the custom Velocity Verlet scheme we use with FLARE, the kinetic energy is one step out of phase with the potential energy. 
To make it easier to analyze the data, we included a few tools to do so (see :class:`Analyze <qmcblip.tools.Analyze>`). 
These tools also align the potential and kinetic energy:

>>> from caf.tools import Analyze
>>> 
>>> # Import the thio.out file
>>> data = Analyze('H2.out')
>>> 
>>> # Create the traj.xyz file
>>> data.to_xyz()
>>> 
>>> # Returns a dictionary containing the keys 'times', 'potential energy', 'kinetic energy',
>>> # 'total energy' and 'temperature'
>>> results = data.get_data()
>>> 
>>> # Plot the energy and save it to energy.png
>>> data.plot_energy(filename="energy.png")

This analyzing tool is only suitable for simulations that were performed with FLARE, and not pure QMC simulations.