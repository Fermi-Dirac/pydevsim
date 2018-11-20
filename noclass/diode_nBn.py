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
import pydevsim
from ds import *
from simdir.physics.new_physics import *
from simdir.physics.ramp2 import *

import numpy as np
# import matplotlib
from matplotlib import pyplot as plt
import ds

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
####
#### Meshing
####
def createMesh(device, region):
  create_1d_mesh(mesh="dio")
  add_1d_mesh_line(mesh="dio", pos=0, ps=1e-7, tag="top")
  add_1d_mesh_line(mesh="dio", pos=0.5e-5, ps=1e-9, tag="mid")
  add_1d_mesh_line(mesh="dio", pos=1e-5, ps=1e-7, tag="bot")
  add_1d_contact  (mesh="dio", name="top", tag="top", material="metal")
  add_1d_contact  (mesh="dio", name="bot", tag="bot", material="metal")
  add_1d_region   (mesh="dio", material="Si", region=region, tag1="top", tag2="bot")
  finalize_mesh(mesh="dio")
  create_device(mesh="dio", device=device)

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
    plt.show()
###### Start #####
# device="MyDevice"
# region="MyRegion"
#
# createMesh(device, region)
ECE_NAME = "ElectronContinuityEquation"
HCE_NAME = "HoleContinuityEquation"

device = "nBn"
meshname = 'nBn_mesh'
n1_region = "n1_region"
b_region = "Barrier"
n2_region = "n2_region"
regions = [n1_region, b_region, n2_region]

n1_thick = 1E-4 # 1e-5  # 0.1 micron, 100 nm
b_thick = 0.5E-4# 0.5e-5
n2_thick = 10E-4 # 10e-5
spacing = min([n1_thick, b_thick, n2_thick])*0.1
b_doping = 5E15
n_eg = 0.121 # eV
n_affinity = 0.052 # less than chi of InAs
b_eg = 1.636 # eV
b_affinity = 1.569 # less than chi of InAs
affinities = [n_affinity, b_affinity, n_affinity]
n_band_gap = 0.121 # eV
b_band_gap = 1.636 # eV

# Devsim requires you create the entire mesh first include interfaces before specifying regions
# Create mesh

ds.create_1d_mesh(mesh=meshname)
# n1
ds.add_1d_mesh_line(mesh=meshname, pos=0, ps=spacing, tag="top_n1")
# add_1d_mesh_line(mesh=meshname, pos=n1_thick/2, ps=spacing, tag="mid_n1")
ds.add_1d_mesh_line(mesh=meshname, pos=n1_thick, ps=spacing, tag="bot_n1")
ds.add_1d_region(mesh=meshname, material="InAs", region=n1_region, tag1="top_n1", tag2="bot_n1")
ds.add_1d_contact(mesh=meshname, name="top_n1", tag="top_n1", material="metal")
# b
# add_1d_mesh_line(mesh=meshname, pos=n1_thick, ps=spacing, tag="top_b")
ds.add_1d_mesh_line(mesh=meshname, pos=n1_thick+(b_thick/2), ps=spacing, tag="mid_b")
ds.add_1d_mesh_line(mesh=meshname, pos=n1_thick+b_thick, ps=spacing, tag="bot_b")
ds.add_1d_region(mesh=meshname, material="InAs", region=b_region, tag1="bot_n1", tag2="bot_b")
ds.add_1d_interface(mesh=meshname, tag="bot_n1", name=f"{n1_region}_to_{b_region}")
# add_1d_interface(mesh=meshname, tag="bot_n1", name=f"{n1_region}_to_{b_region}")
# CreateInterfaceModel(device=device, interface=)
# n_2
# add_1d_mesh_line(mesh=meshname, pos=n1_thick+b_thick, ps=spacing, tag="top_n2")
ds.add_1d_mesh_line(mesh=meshname, pos=n1_thick+b_thick+(n2_thick/2), ps=spacing, tag="mid_n2")
ds.add_1d_mesh_line(mesh=meshname, pos=(n1_thick+b_thick+n2_thick), ps=spacing, tag="bot_n2")
ds.add_1d_region(mesh=meshname, material="InAs", region=n2_region, tag1="bot_b", tag2="bot_n2")

ds.add_1d_interface(mesh=meshname, tag="bot_b", name=f"{b_region}_to_{n2_region}")
ds.add_1d_contact(mesh=meshname, name="bot_n2", tag="bot_n2", material="metal")
ds.finalize_mesh(mesh=meshname)
ds.create_device(mesh=meshname, device=device)

####
#### Set parameters for 300 K
####
ds.set_parameter(name="T", value=300)
for region in regions:
    if region == b_region:
        affinity = b_affinity
        bandgap=b_band_gap
        doping = b_doping
    else:
        affinity = n_affinity
        bandgap = n_band_gap
        doping = 1E11
    SetSiliconParameters(device, region, Affinity=affinity, EG300=bandgap)
    set_parameter(device=device, region=region, name="taun", value=1e-8)
    set_parameter(device=device, region=region, name="taup", value=1e-8)
    set_parameter(device=device, region=region, name="n1", value=1e10)
    set_parameter(device=device, region=region, name="p1", value=1e10)


    ####
    #### NetDoping
    ####
    CreateNodeModel(device, region, "Acceptors", f"{1e11:.6e}")# "1.0e18*step(0.5e-5-x)")
    CreateNodeModel(device, region, "Donors",    f"{doping:.6e}")#"1.0e18*step(x-0.5e-5)")
    CreateNodeModel(device, region, "NetDoping", "Donors-Acceptors")
    # print_node_values(device=device, region=region, name="NetDoping")

    ####
    #### Create Potential, Potential@n0, Potential@n1
    ####
    CreateSolution(device, region, "Potential")

    ####
    #### Create potential only physical models
    ####
    CreateSiliconPotentialOnly(device, region)
    for contact in ds.get_contact_list(device=device):
        CreateSiliconPotentialOnlyContact(device, region, contact)

    ####
    #### Set up the contacts applying a bias
    ####
    for i in get_contact_list(device=device):
      CreateSiliconPotentialOnlyContact(device, region, i)
for i in get_contact_list(device=device):
  set_parameter(device=device, name=GetContactBiasName(i), value=0.0)
    # ####
    # #### Initial DC solution
    # ####
solve(type="dc", absolute_error=1.0, relative_error=1e-10, maximum_iterations=30)
for region in regions:
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
pydevsim.get_ds_status()
pydevsim.plot_charge(regions=regions)
solve(type="dc", solver_type='iterative', absolute_error=1e1, relative_error=1e-10, maximum_iterations=500, info=True )
print(">>After Equil<<")
pydevsim.get_ds_status()
pydevsim.plot_charge(regions=regions)
####
#### Ramp the bias to 0.5 Volts
####
volts = []
current = []
# v = -0.5
for v in np.linspace(-1, 1, 10):
  v = float(v)
  set_parameter(device=device, name=GetContactBiasName("top_n1"), value=v)
  solve(type="dc", absolute_error=1e1, relative_error=1e-2, maximum_iterations=200)
  # PrintCurrents(device, "top_n1")
  # PrintCurrents(device, "bot_n2")
  e_current =  abs(ds.get_contact_current(device=device, contact='top_n1', equation=ECE_NAME))
  e_current += abs(ds.get_contact_current(device=device, contact='top_n1', equation=HCE_NAME))
  e_current += abs(ds.get_contact_current(device=device, contact='bot_n2', equation=ECE_NAME))
  e_current += abs(ds.get_contact_current(device=device, contact='bot_n2', equation=HCE_NAME))
  current.append(e_current)
  volts.append(v)
  pydevsim.plot_potential(regions=regions)
  # v += 0.1
plt.figure()
plt.plot(volts, current)
  # pydevsim.plot_charge(regions=regions)
  # pydevsim.plot_current(regions=regions)
plt.show()

ds.write_devices(file="diode_1d.tec", type="tecplot")
plt.clf()
# edge_average_model(device=device, region=region, node_model="x", edge_model="xmid")
# xmid=get_edge_model_values(device=device, region=region, name="xmid")
#efields = ("Jn_arora_lf", "Jp_arora_lf" )
#efields = ("Jn", "Jp", "Jn_arora_lf", "Jp_arora_lf" )

efields = ("Jn", "Jp")
pydevsim.plot_current(regions=regions, current_names=efields)
# y=get_edge_model_values(device=device, region=region, name=efields[0])
# ymin=min(y)
# ymax=max(y)
# for i in efields:
#   y=get_edge_model_values(device=device, region=region, name=i)
#   if min(y) < ymin:
#     ymin = min(y)
#   elif max(y) > ymax:
#     ymax = max(y)
#   plt.plot(xmid, y)
# plt.xlabel('x (cm)')
# plt.ylabel('J (A/cm^2)')
# plt.legend(efields)
# # plt.axis([min(x), max(x), 0.5*ymin, 2*ymax])
# plt.savefig("diode_1d_current.png")
# plt.show()
# print (ymin)
# print (ymax)

plt.clf()
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
  plt.plot(xmid, y)
plt.xlabel('x (cm)')
plt.ylabel('J (A/cm^2)')
plt.legend(efields)
plt.axis([min(x), max(x), 0.5*ymin, 2*ymax])
plt.savefig("diode_1d_mobility.png")
plt.show()
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
    plt.semilogy(np.array(x)*10000, y)
plt.xlabel('x (nm)')
plt.ylabel('Density (#/cm^3)')
plt.legend(fields)
ymax *= 10
plt.axis([min(x)*10000, max(x)*10000, ymin, ymax])
plt.savefig("USRH.png")
plt.show()

