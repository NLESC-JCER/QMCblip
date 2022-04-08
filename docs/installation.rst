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


Installing CHAMP
=================

To use this software, you need a special version of CHAMP which can export the forces and energy to a file. 
For this you need to ase-coupling_ branch. 
Simply clone this and build in the usual way, as specified by their documentation.

.. _ase-coupling: https://github.com/filippi-claudia/champ/tree/ase-coupling
