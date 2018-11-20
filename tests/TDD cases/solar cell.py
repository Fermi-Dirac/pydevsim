from matplotlib import pyplot as plt
from pydevsim import materials, device, environment

# TDD case for a resistor made of silicon given a particular doping profile

p_Si = materials.Silicon(thickness=10, nd=1E15)  # Si constructor handles everything we want by default. required parameter is thickness
n_Si = materials.Silicon(thickness=10, na=1E15)  # Si constructor handles everything we want by default. required parameter is thickness
sun = environment.Sunlight(spectra='AM1.5', power=1000)
outdoors = environment(light=sun)  # default constructor, V=0, T=293 K, B=0, E=0, adding light source
solar_cell = device.Device(material_list=[p_Si, n_Si], interfaces=None, contacts=None, R_series=1, environment=outdoors)
volts, curr = range(-5,5,1), []
measurement = solar_cell.sweep_iv()
# this method sweeps the IV and saves salient model paraemters and data into a measurement object which contains
# pandas dataframes

# or do it manually:
for volt in volts:
    outdoors.voltage = volt
    solar_cell.solve()
    curr.append(solar_cell.current)

plt.plot(volts, curr)
plt.show()