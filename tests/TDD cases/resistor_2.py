from matplotlib import pyplot as plt
from pydevsim import material, device, region

# TDD case for a resistor made of silicon given a particular doping profile
sim_space = Simulation(meshname='lol', device_name='lol')
material.Silicon(nd=1E15)
region.Region1D(sim=sim_space, material=p_Si)
device.Device1D(sim=simspace, region)

resistor = device.Device1D(name='resistor', )
resistor.add_region()

Si = material.Silicon(Nd=1E15)  # Si constructor handles everything we want by default. required parameter is thickness
p_Si = region.Region1D(material=Si, thickness=1)
resistor = device.Device1D(regions=[p_Si], R_series=1, )
volts, curr = range(-5,5,1), []
measurement = resistor.sweep_iv()
# this method sweeps the IV and saves salient model paraemters and data into a measurement object which contains
# pandas dataframes

# or do it manually:
for volt in volts:
    resistor.set=volt
    resistor.solve()
    curr.append(resistor.total_current)

plt.plot(volts, curr)
plt.show()