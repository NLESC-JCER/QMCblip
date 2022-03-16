import numpy as np
from flare import otf_parser
import matplotlib.pyplot as plt


class Analyze():

    def __init__(self, file):
        self.otf = otf_parser.OtfAnalysis(file, calculate_energy=True)
        self.results = None

    def get_data(self):
        dt = self.otf.header['dt']

        frames = np.arange(len(self.otf.times)+1)

        self.results = dict()

        self.dft = self.otf.dft_frames
        self.nondft = [i for i in frames if i not in self.dft]

        totE = np.zeros(len(frames))
        kinE = np.zeros(len(frames))
        potE = np.zeros(len(frames))
        temp = np.zeros(len(frames))

        dftframe = 0
        for i in range(len(frames)):
            if i in self.dft:
                totE[i] = self.otf.gp_thermostat['total energy'][dftframe]
                kinE[i] = self.otf.gp_thermostat['kinetic energy'][dftframe]
                potE[i] = self.otf.gp_thermostat['potential energy'][dftframe]
                temp[i] = self.otf.gp_thermostat['temperature'][dftframe]
                dftframe += 1
            else:
                totE[i] = self.otf.thermostat['total energy'][i-1]
                kinE[i] = self.otf.thermostat['kinetic energy'][i-1]
                potE[i] = self.otf.thermostat['potential energy'][i-1]
                temp[i] = self.otf.thermostat['temperature'][i-1]


        self.results['times'] = frames[0:-2] * dt
        self.results['potential energy'] = potE[1:-1]
        self.results['kinetic energy'] = kinE[2:]
        self.results['total energy'] = kinE[2:] + potE[1:-1]
        self.results['temperature'] = temp[2:]

        if (max(self.dft) > max(self.nondft)):
            self.dft = self.dft[:-1]
        else:
            self.nondft = self.nondft[:-1]

        return self.results
    
    def plot_energy(self, filename=""):
        if self.results is None:
            self.get_data()

        times = self.results['times']
        totE = self.results['total energy']
        potE = self.results['potential energy']
        kinE = self.results['kinetic energy']
        plt.figure()
        plt.grid()
        plt.title("Energy")
        plt.plot(times, totE - totE[0], label="Total Energy")
        plt.plot(times, potE - potE[0], label="Potential Energy")
        plt.plot(times, kinE - kinE[0], label="Kinetic Energy")
        plt.scatter([times[index-1] for index in self.dft[1:]], [totE[index-1] - totE[0] for index in self.dft[1:]], marker='x', label="DFT Calculations", color='black')
        #plt.scatter([times[index-1] for index in self.nondft], [totE[index-1] - totE[0] for index in self.nondft], marker='x', label="FLARE Calculations", color='blue')
        plt.legend()
        plt.xlabel("Time (ps)")
        plt.ylabel("Energy (eV)")
        if (filename != ""):
            plt.savefig(filename)

    def plot_temperature(self, filename=""):
        if self.results is None:
            self.get_data()

            plt.figure()

        times = self.results['times']
        temp = self.results['temperature']
        plt.figure()
        plt.grid()
        plt.title("Temperature")
        plt.plot(times, temp)
        plt.xlabel("Time (ps)")
        plt.ylabel("Temperature (K)")
        if (filename != ""):
            plt.savefig(filename)