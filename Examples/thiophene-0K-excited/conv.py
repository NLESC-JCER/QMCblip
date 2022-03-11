
import numpy as np
from flare import otf_parser
import matplotlib.pyplot as plt


otf_object = otf_parser.OtfAnalysis('thio.out',  calculate_energy=True)
otf_object.to_xyz('traj.xyz')

totE = otf_object.thermostat['total energy']
kinE = otf_object.thermostat['kinetic energy']
potE = otf_object.thermostat['potential energy']
times = 0.5e-3 * np.array([i for i in range(len(totE))])

totE = (np.array(totE) - totE[0])
potE = (np.array(potE) - potE[0])
kinE = (np.array(kinE) - kinE[0])

dft = otf_object.dft_frames
nondft = [i for i in [j for j in range(len(times))] if i not in dft]

plt.figure()
plt.grid()
plt.title("Energy")
plt.plot(times, totE, label="Total Energy")
plt.plot(times, potE, label="Potential Energy")
plt.plot(times, kinE, label="Kinetic Energy")
plt.scatter([times[index] for index in dft], [totE[index] for index in dft], label="DFT Calculations", color='black')
#plt.scatter([times[index] for index in nondft], [totE[index] for index in nondft], label="FLARE Calculations", color='blue')
plt.legend()
plt.xlabel("Time (ps)")
plt.ylabel("Energy (eV)")
plt.savefig("energy.png")
plt.show()