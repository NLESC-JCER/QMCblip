Installation
=========================

Via Python Package
-----------------------------

The latest release of QMCblip can be installed usint the pypi package manager with:

``pip install qmcblip`` 


Via GitHub
-------------

For user who would like to contribute, the code is hosted on GitHub_.

.. _GitHub: https://github.com/NLESC-JCER/QMCblip

To install the code

 * clone the repository ``git clone https://github.com/NLESC-JCER/QMCblip.git``
 * go there ``cd QMCblip``
 * install the module ``pip install -e ./``


Installing Dependencies
=======================

Installing CHAMP
----------------

To use this software, you need a special version of CHAMP which can export the forces and energy to a file. 
For this you need to ase-coupling_ branch. 
Simply clone this and build in the usual way, as specified by their documentation.

.. _ase-coupling: https://github.com/filippi-claudia/champ/tree/ase-coupling

Installing GAMESS
-----------------

GAMESS can be found at gamess_. To make CHAMP and GAMESS work together, we need to change some settings.
First, make sure that GAMESS is on your path, specifically the rungms file.
Secondly, in the ``gms-files.csh`` file change the ``setenv EXTBAS /dev/null`` to ``setenv EXTBAS /your/champ/location/pool/BFD/BASIS_gamess/BFD_Basis.EXTBAS``. Now GAMESS can use the BFD basis.

.. _gamess: https://www.msg.chem.iastate.edu/gamess/

