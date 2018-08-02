from pydevsim import material, device, environment, mesh
from matplotlib import pyplot as plt
# TDD case for a resistor made of silicon given a particular doping profile

Si = material.Silicon(thickness=10, nd=1E15)  # Si constructor handles everything we want by default. required parameter is thickness
stp = environment.Environment()  # default constructor, V=0, T=293 K, B=0, E=0, no light source
resistor = device.Device1D(material_list=[Si], interfaces=None, contacts=None, R_series=1, environment=stp)
volts, curr = range(-5,5,1), []
measurement = resistor.sweep_iv()
# this method sweeps the IV and saves salient model paraemters and data into a measurement object which contains
# pandas dataframes

# or do it manually:
for volt in volts:
    resistor.voltage=volt
    resistor.solve()
    curr.append(resistor.current)

plt.plot(volts, curr)
plt.show()