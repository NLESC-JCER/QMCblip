Analyzing results
-----------------

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