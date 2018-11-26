# Copyright 2016 Devsim LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from ds import *
from simdir.physics.new_physics import *
from simdir.physics.ramp2 import *

import numpy as np
# import matplotlib
import matplotlib.pyplot
from matplotlib import pyplot as plt
import ds
ECE_NAME = "ElectronContinuityEquation"
HCE_NAME = "HoleContinuityEquation"
#####
# dio1
#
# Make doping a step function
# print dat to text file for viewing in grace
# verify currents analytically
# in dio2 add recombination
#
def get_ds_status():
    devices = ds.get_device_list()
    for device in devices:
        print("Device: " + device)
        regions = ds.get_region_list(device=device)
        for region in regions:
            print("\tRegion :" + region)
            params = ds.get_parameter_list(device=device, region=region)
            for param in params:
                val = ds.get_parameter(device=device, region=region, name=param)
                print(f"\t\t{param} = {val}")
            n_models = ds.get_node_model_list(device=device, region=region)
            for node_model in n_models:
                nmvals = ds.get_node_model_values(device=device, region=region, name=node_model)
                print(f"\t\t Node Model '{node_model}' = {nmvals!s}")
            e_models = ds.get_edge_model_list(device=device, region=region)
            for edge_model in e_models:
                emvals = ds.get_edge_model_values(device=device, region=region, name=edge_model)
                print(f"\t\t Edge Model '{edge_model}' = {emvals!s}")
        contacts = ds.get_contact_list(device=device)
        for contact in contacts:
            print("\tContact : " + contact)
            c_eqs = ds.get_contact_equation_list(device=device, contact=contact)
            for ceq in c_eqs:
                print("\t\tContact Equation : " + ceq)


def plot_charge(regions):
    fields = ("Electrons", "Holes", "Donors", "Acceptors")
    plt.figure()
    for var_name in fields:
        total_x = np.array([])
        total_y = []
        for region in regions:
            x = np.array(get_node_model_values(device=device, region=region, name="x"))
            # print(var_name, min(x), max(x), min(x)*1e4, max(x)*1e4)
            total_x = np.append(total_x, x)
            y=get_node_model_values(device=device, region=region, name=var_name)
            total_y.extend(y)
        # plt.axis([min(x), max(x), ymin, ymax])
        plt.semilogy(np.array(total_x)*1e4, total_y)
    plt.xlabel('x (um)')
    plt.ylabel('Density (#/cm^3)')
    plt.legend(fields)
    # plt.savefig("diode_1d_density.png")
####
#### Meshing
####
def createMesh(device, region):
  create_1d_mesh(mesh="dio")
  add_1d_mesh_line(mesh="dio", pos=0, ps=5e-8, tag="top")
  add_1d_mesh_line(mesh="dio", pos=0.5e-5, ps=1e-9, tag="mid")
  add_1d_mesh_line(mesh="dio", pos=1e-5, ps=5e-8, tag="bot")
  add_1d_contact  (mesh="dio", name="top", tag="top", material="metal")
  add_1d_contact  (mesh="dio", name="bot", tag="bot", material="metal")
  add_1d_region   (mesh="dio", material="Si", region=region, tag1="top", tag2="bot")
  finalize_mesh(mesh="dio")
  create_device(mesh="dio", device=device)

device="MyDevice"
region="MyRegion"

createMesh(device, region)

####
#### Set parameters for 300 K
####
set_parameter(name="T", value=300)
SetSiliconParameters(device, region)
set_parameter(device=device, region=region, name="taun", value=1e-8)
set_parameter(device=device, region=region, name="taup", value=1e-8)
set_parameter(device=device, region=region, name="n1", value=1e10)
set_parameter(device=device, region=region, name="p1", value=1e10)


####
#### NetDoping
####
CreateNodeModel(device, region, "Acceptors", "1.0e18*step(0.5e-5-x)")
CreateNodeModel(device, region, "Donors",    "1.0e18*step(x-0.5e-5)")
CreateNodeModel(device, region, "NetDoping", "Donors-Acceptors")
print_node_values(device=device, region=region, name="NetDoping")

####
#### Create Potential, Potential@n0, Potential@n1
####
CreateSolution(device, region, "Potential")

####
#### Create potential only physical models
####
CreateSiliconPotentialOnly(device, region)

####
#### Set up the contacts applying a bias
####
for i in get_contact_list(device=device):
  set_parameter(device=device, name=GetContactBiasName(i), value=0.0)
  CreateSiliconPotentialOnlyContact(device, region, i)


####
#### Initial DC solution
####
solve(type="dc", absolute_error=1.0e2, relative_error=1e-1, maximum_iterations=100, solver_type='iterative')

####
#### drift diffusion solution variables
####
CreateSolution(device, region, "Electrons")
CreateSolution(device, region, "Holes")

CreateEField(device, region)
CreateDField(device, region)
opts = CreateAroraMobilityLF(device, region)
opts = CreateHFMobility(device, region, **opts)
# CreateHFMobility(device, region, **opts)

set_parameter(device=device, region=region, name="BETAN",  value=2.0)
set_parameter(device=device, region=region, name="BETAP",  value=1.0)
set_parameter(device=device, region=region, name="VSATN0",  value=2.4e7)
set_parameter(device=device, region=region, name="VSATP0",  value=2.4e7)
set_parameter(device=device, region=region, name="VSATN.A",  value=0.8)
set_parameter(device=device, region=region, name="VSATP.A",  value=0.8)


####
#### create initial guess from dc only solution
####
set_node_values(device=device, region=region, name="Electrons", init_from="IntrinsicElectrons")
set_node_values(device=device, region=region, name="Holes",     init_from="IntrinsicHoles")

# import physics.model_create
#physics.model_create.debug=True
###
### Set up equations
###
CreateSiliconDriftDiffusion(device, region, **opts)
for i in get_contact_list(device=device):
  CreateSiliconDriftDiffusionContact(device, region, i, Jn=opts['Jn'], Jp=opts['Jp'])

###
### Drift diffusion simulation at equilibrium
###
# set_parameter(device=device, name=GetContactBiasName("bot"), value=2.0)
solve(type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=100)
print(">>After Equil<<")

get_ds_status()
####
#### Ramp the bias to 0.5 Volts
####
volts=[]
current=[]
v = 0.0
while v < 0.51:
  # plot_charge([region])
  set_parameter(device=device, name=GetContactBiasName("bot"), value=v)
  set_parameter(device=device, name=GetContactBiasName("top"), value=0.0)
  print(v)
  solve(type="dc", absolute_error=1e6, relative_error=1e-12, maximum_iterations=300)
  # PrintCurrents(device, "top")
  # PrintCurrents(device, "bot")
  for contact in ds.get_contact_list(device=device):
    e_current = abs(ds.get_contact_current(device=device, contact=contact, equation=ECE_NAME))
    e_current += abs(ds.get_contact_current(device=device, contact=contact, equation=HCE_NAME))
  current.append(e_current)
  volts.append(v)
  v += 0.02
# plt.show()
plt.figure()
plt.semilogy(volts, current)
plt.figure()
plt.plot(volts, current)
  # pydevsim.plot_charge(regions=regions)
  # pydevsim.plot_current(regions=regions)
plt.show()

write_devices(file="diode_1d.tec", type="tecplot")

x=get_node_model_values(device=device, region=region, name="x")
ymax = 10
ymin = 10
fields = ("Electrons", "Holes", "Donors", "Acceptors")
for i in fields:
    y=get_node_model_values(device=device, region=region, name=i)
    if (max(y) > ymax):
      ymax = max(y)
    matplotlib.pyplot.semilogy(x, y)
matplotlib.pyplot.xlabel('x (cm)')
matplotlib.pyplot.ylabel('Density (#/cm^3)')
matplotlib.pyplot.legend(fields)
ymax *= 10
matplotlib.pyplot.axis([min(x), max(x), ymin, ymax])
matplotlib.pyplot.savefig("diode_1d_density.png")
matplotlib.pyplot.show()

matplotlib.pyplot.clf()
edge_average_model(device=device, region=region, node_model="x", edge_model="xmid")
xmid=get_edge_model_values(device=device, region=region, name="xmid")
#efields = ("Jn_arora_lf", "Jp_arora_lf" )
#efields = ("Jn", "Jp", "Jn_arora_lf", "Jp_arora_lf" )
efields = ("Jn", "Jp")
y=get_edge_model_values(device=device, region=region, name=efields[0])
ymin=min(y)
ymax=max(y)
for i in efields:
  y=get_edge_model_values(device=device, region=region, name=i)
  if min(y) < ymin:
    ymin = min(y)
  elif max(y) > ymax:
    ymax = max(y)
  matplotlib.pyplot.plot(xmid, y)
matplotlib.pyplot.xlabel('x (cm)')
matplotlib.pyplot.ylabel('J (A/cm^2)')
matplotlib.pyplot.legend(efields)
matplotlib.pyplot.axis([min(x), max(x), 0.5*ymin, 2*ymax])
matplotlib.pyplot.savefig("diode_1d_current.png")
matplotlib.pyplot.show()
print (ymin)
print (ymax)

matplotlib.pyplot.clf()
edge_average_model(device=device, region=region, node_model="x", edge_model="xmid")
xmid=get_edge_model_values(device=device, region=region, name="xmid")
efields = ("mu_arora_n_lf", "mu_arora_p_lf", "mu_n", "mu_p",  )
#efields = ("Jn", "Jp", "Jn_arora_lf", "Jp_arora_lf" )
y=get_edge_model_values(device=device, region=region, name=efields[0])
ymin=min(y)
ymax=max(y)
for i in efields:
  y=get_edge_model_values(device=device, region=region, name=i)
  if min(y) < ymin:
    ymin = min(y)
  elif max(y) > ymax:
    ymax = max(y)
  matplotlib.pyplot.plot(xmid, y)
matplotlib.pyplot.xlabel('x (cm)')
matplotlib.pyplot.ylabel('J (A/cm^2)')
matplotlib.pyplot.legend(efields)
matplotlib.pyplot.axis([min(x), max(x), 0.5*ymin, 2*ymax])
matplotlib.pyplot.savefig("diode_1d_mobility.png")
matplotlib.pyplot.show()
print (ymin)
print (ymax)


#x=get_node_model_values(device=device, region=region, name="x")
ymax = 10
ymin = 10
fields = ("USRH",)
for i in fields:
    y=get_node_model_values(device=device, region=region, name=i)
    if (max(y) > ymax):
      ymax = max(y)
    matplotlib.pyplot.semilogy(np.array(x)*10000, y)
matplotlib.pyplot.xlabel('x (nm)')
matplotlib.pyplot.ylabel('Density (#/cm^3)')
matplotlib.pyplot.legend(fields)
ymax *= 10
matplotlib.pyplot.axis([min(x)*10000, max(x)*10000, ymin, ymax])
matplotlib.pyplot.savefig("USRH.png")
matplotlib.pyplot.show()
