from caf.tools import Analyze
import matplotlib.pyplot as plt

t = Analyze('thio.out')

t.to_xyz()

t.get_data()

t.plot_energy("energy.png")

t.plot_temperature('temperature.png')