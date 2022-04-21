Creating Wavefunction
---------------------

Important warning: this part is currently still experimental and only confirmed to work for C2.

QMCblip also includes some tools to make the wavefunction for CHAMP using GAMESS. To use it, first import:

>>> from qmcblip.gamess.utils import WavefunctionCreator

Also create your molecular system:

>>> from ase import Atoms
>>> atoms = Atoms('C2', [(0,0,-0.7), (0,0, 0.7)])

We can now initialize our :obj:`WavefunctionCreator <qmcblip.gamess.utils.WavefunctionCreator>` object. The first argument are our atoms and the second argument is the location of the CHAMP base directory. This is not the directory containing ``vmc.mov1``, but the directory one level up.

>>> wf = WavefunctionCreator(atoms, '../../champ')

We can now run the GAMESS simulations. See the `GAMESS documentation`_ for more information. Using keywords and dictionaries you can supply additional GAMESS settings or change the defaults. Using the userscr keyword you can supply the location of the USERSCR directory. Usually the ASE GAMESS calculator can find it on its own, but if not you need to supply it.

>>> wf.setup_rhf()
>>> wf.setup_cas(system=dict(mwords=500), drt=dict(nmcc=2, ndoc=2, nval=2))
>>> wf.setup_ci(system=dict(mwords=500), cidrt=dict(nfzc=2, ndoc=2, nval=2))

After GAMESS is done we can convert the GAMESS output to CHAMP wavefunction.

>>> wf.convert_to_champ()

And we can also quickly make a CHAMP input file (``vmc.inp``).

>>> input = wf.create_champ_input()

Here ``input`` is of the :obj:`Settings <qmcblip.champio.Settings>` type. With that we can easily run a CHAMP simulation:

>>> from qmcblip.champ import CHAMP
>>> atoms.calc = CHAMP(champ_loc='../../champ/bin/vmc.mov1', settings = input)

For more examples, see :doc:`../examples/examples`.

.. _`GAMESS documentation`: https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/docs-input.txt