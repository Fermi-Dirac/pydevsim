# from ds import *
from ds import create_1d_mesh, add_1d_mesh_line, add_1d_region, add_1d_contact, finalize_mesh, create_device
# from python_packages.simple_physics import *
from ds import print_node_values, solve, get_node_model_values, edge_average_model, get_edge_model_values, get_contact_current
from python_packages.simple_physics import SetSiliconParameters, CreateNodeModel, CreateSolution, \
    CreateSiliconPotentialOnly, get_contact_list, CreateSiliconPotentialOnlyContact, set_parameter, GetContactBiasName
from python_packages.simple_physics import PrintCurrents, write_devices, set_node_values, CreateSiliconDriftDiffusionAtContact, CreateSiliconDriftDiffusion
from python_packages.simple_physics import ece_name, hce_name
from matplotlib import pyplot as plt
import numpy as np

def CreateMesh(device, region, meshname='diode'):
    '''
      Meshing
    '''
    create_1d_mesh(mesh=meshname)
    add_1d_mesh_line(mesh=meshname, pos=0, ps=1e-7, tag="top")
    add_1d_mesh_line(mesh=meshname, pos=0.5e-5, ps=1e-9, tag="mid")
    add_1d_mesh_line(mesh=meshname, pos=1e-5, ps=1e-7, tag="bot")
    add_1d_contact(mesh=meshname, name="top", tag="top", material="metal")
    add_1d_contact(mesh=meshname, name="bot", tag="bot", material="metal")
    add_1d_region(mesh=meshname, material="Si", region=region, tag1="top", tag2="bot")
    finalize_mesh(mesh=meshname)
    create_device(mesh=meshname, device=device)


# this is the mesh for the ssac device
# TODO: use CreateMesh and update regressions
def CreateMesh2(device, region):
    create_1d_mesh(mesh="dio")
    add_1d_mesh_line(mesh="dio", pos=0, ps=1e-7, tag="top")
    add_1d_mesh_line(mesh="dio", pos=0.5e-5, ps=1e-8, tag="mid")
    add_1d_mesh_line(mesh="dio", pos=1e-5, ps=1e-7, tag="bot")
    add_1d_contact(mesh="dio", name="top", tag="top", material="metal")
    add_1d_contact(mesh="dio", name="bot", tag="bot", material="metal")
    add_1d_region(mesh="dio", material="Si", region=region, tag1="top", tag2="bot")
    finalize_mesh(mesh="dio")
    create_device(mesh="dio", device=device)


# def Create2DMesh(device, region):
#   create_2d_mesh(mesh="dio")
#   add_2d_mesh_line(mesh="dio", dir="x", pos=0,      ps=1e-6)
#   add_2d_mesh_line(mesh="dio", dir="x", pos=0.5e-5, ps=1e-8)
#   add_2d_mesh_line(mesh="dio", dir="x", pos=1e-5,   ps=1e-6)
#   add_2d_mesh_line(mesh="dio", dir="y", pos=0,      ps=1e-6)
#   add_2d_mesh_line(mesh="dio", dir="y", pos=1e-5,   ps=1e-6)
#
#   add_2d_mesh_line(mesh="dio", dir="x", pos=-1e-8,    ps=1e-8)
#   add_2d_mesh_line(mesh="dio", dir="x", pos=1.001e-5, ps=1e-8)
#
#   add_2d_region(mesh="dio", material="Si", region=region)
#   add_2d_region(mesh="dio", material="Si", region="air1", xl=-1e-8,  xh=0)
#   add_2d_region(mesh="dio", material="Si", region="air2", xl=1.0e-5, xh=1.001e-5)
#
#   add_2d_contact(mesh="dio", name="top", material="metal", region=region, yl=0.8e-5, yh=1e-5, xl=0, xh=0, bloat=1e-10)
#   add_2d_contact(mesh="dio", name="bot", material="metal", region=region, xl=1e-5,   xh=1e-5, bloat=1e-10)
#
#   finalize_mesh(mesh="dio")
#   create_device(mesh="dio", device=device)

# def Create2DGmshMesh(device, region):
#   #this reads in the gmsh format
#   create_gmsh_mesh (mesh="diode2d", file="gmsh_diode2d.msh")
#   add_gmsh_region  (mesh="diode2d", gmsh_name="Bulk",    region=region, material="Silicon")
#   add_gmsh_contact (mesh="diode2d", gmsh_name="Base",    region=region, material="metal", name="top")
#   add_gmsh_contact (mesh="diode2d", gmsh_name="Emitter", region=region, material="metal", name="bot")
#   finalize_mesh    (mesh="diode2d")
#   create_device    (mesh="diode2d", device=device)
#
# def Create3DGmshMesh(device, region):
#   #this reads in the gmsh format
#   create_gmsh_mesh (mesh="diode3d", file="gmsh_diode3d.msh")
#   add_gmsh_region  (mesh="diode3d", gmsh_name="Bulk",    region=region, material="Silicon")
#   add_gmsh_contact (mesh="diode3d", gmsh_name="Base",    region=region, material="metal", name="top")
#   add_gmsh_contact (mesh="diode3d", gmsh_name="Emitter", region=region, material="metal", name="bot")
#   finalize_mesh    (mesh="diode3d")
#   create_device    (mesh="diode3d", device=device)


def SetParameters(device, region):
    '''
      Set parameters for 300 K
    '''
    SetSiliconParameters(device, region, 300)


def SetNetDoping(device, region):
    '''
      NetDoping
    '''
    CreateNodeModel(device, region, "Acceptors", "1.0e18*step(0.5e-5-x)")
    CreateNodeModel(device, region, "Donors", "1.0e18*step(x-0.5e-5)")
    CreateNodeModel(device, region, "NetDoping", "Donors-Acceptors")


def InitialSolution(device, region, circuit_contacts=None):
    # Create Potential, Potential@n0, Potential@n1
    CreateSolution(device, region, "Potential")

    # Create potential only physical models
    CreateSiliconPotentialOnly(device, region)

    # Set up the contacts applying a bias
    for i in get_contact_list(device=device):
        if circuit_contacts and i in circuit_contacts:
            CreateSiliconPotentialOnlyContact(device, region, i, True)
        else:
            ###print "FIX THIS"
            ### it is more correct for the bias to be 0, and it looks like there is side effects
            set_parameter(device=device, name=GetContactBiasName(i), value=0.0)
            CreateSiliconPotentialOnlyContact(device, region, i)


def DriftDiffusionInitialSolution(device, region, circuit_contacts=None):
    ####
    #### drift diffusion solution variables
    ####
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")

    ####
    #### create initial guess from dc only solution
    ####
    set_node_values(device=device, region=region, name="Electrons", init_from="IntrinsicElectrons")
    set_node_values(device=device, region=region, name="Holes", init_from="IntrinsicHoles")

    ###
    ### Set up equations
    ###
    CreateSiliconDriftDiffusion(device, region)
    for i in get_contact_list(device=device):
        if circuit_contacts and i in circuit_contacts:
            CreateSiliconDriftDiffusionAtContact(device, region, i, True)
        else:
            CreateSiliconDriftDiffusionAtContact(device, region, i)


#####
##### Ramp the bias to 0.5 Volts
#####
# v = 0.0
# while v < 0.51:
#  set_parameter(device=device, name=GetContactBiasName("top"), value=v)
#  solve(type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=30)
#  PrintCurrents(device, "top")
#  PrintCurrents(device, "bot")
#  v += 0.1
#
# write_devices(file="diode_1d.dat", type="tecplot")
##import matplotlib
##import matplotlib.pyplot
##x=get_node_model_values(device=device, region=region, name="x")
##ymax = 10
##ymin = 10
##fields = ("Electrons", "Holes", "Donors", "Acceptors")
##for i in fields:
##    y=get_node_model_values(device=device, region=region, name=i)
##    if (max(y) > ymax):
##      ymax = max(y)
##    matplotlib.pyplot.semilogy(x, y)
##matplotlib.pyplot.xlabel('x (cm)')
##matplotlib.pyplot.ylabel('Density (#/cm^3)')
##matplotlib.pyplot.legend(fields)
##ymax *= 10
##matplotlib.pyplot.axis([min(x), max(x), ymin, ymax])
##matplotlib.pyplot.savefig("diode_1d_density.eps")
##
##matplotlib.pyplot.clf()
##edge_average_model(device=device, region=region, node_model="x", edge_model="xmid")
##xmid=get_edge_model_values(device=device, region=region, name="xmid")
##efields = ("ElectronCurrent", "HoleCurrent", )
##y=get_edge_model_values(device=device, region=region, name="ElectronCurrent")
##ymin=min(y)
##ymax=max(y)
##for i in efields:
##  y=get_edge_model_values(device=device, region=region, name=i)
##  if min(y) < ymin:
##    ymin = min(y)
##  elif max(y) > ymax:
##    ymax = max(y)
##  matplotlib.pyplot.plot(xmid, y)
##matplotlib.pyplot.xlabel('x (cm)')
##matplotlib.pyplot.ylabel('J (A/cm^2)')
##matplotlib.pyplot.legend(efields)
##matplotlib.pyplot.axis([min(x), max(x), 0.5*ymin, 2*ymax])
##matplotlib.pyplot.savefig("diode_1d_current.eps")
##print ymin
##print ymax


device = "MyDevice"
region = "MyRegion"

CreateMesh(device=device, region=region)

SetParameters(device=device, region=region)
set_parameter(device=device, region=region, name="taun", value=1e-8)
set_parameter(device=device, region=region, name="taup", value=1e-8)

SetNetDoping(device=device, region=region)

print_node_values(device=device, region=region, name="NetDoping")

InitialSolution(device, region)

# Initial DC solution
solve(type="dc", absolute_error=1.0, relative_error=1e-10, maximum_iterations=30)

DriftDiffusionInitialSolution(device, region)
###
### Drift diffusion simulation at equilibrium
###
solve(type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=30)

####
#### Ramp the bias to 0.5 Volts
####
# bias_volt = 0.0
# while bias_volt < 0.51:

def get_total_current(device):
    tot_current = {}
    for contact in ['top', 'bot']:
        e_current = get_contact_current(device=device, contact=contact, equation=ece_name)
        h_current = get_contact_current(device=device, contact=contact, equation=hce_name)
        tot_current[contact] = e_current + h_current
    return tot_current['top'], tot_current['bot']

volt_sweep = np.arange(-0.5, .7, 0.01)
top_currents = []
bot_currents = []
for bias_volt in volt_sweep:
# v = 0.0
# while v < 0.51:
#     bias_volt = v
    bias_volt = float(bias_volt)
    # print(f"Bias voltage is {bias_volt} and is of type {type(bias_volt)}")
    set_parameter(device=device, name="{0}_bias".format("top"), value=bias_volt)
    solve(type="dc", absolute_error=1e9, relative_error=1e-10, maximum_iterations=50)
    # PrintCurrents(device, "top")
    # PrintCurrents(device, "bot")
    tot_top, tot_bot = get_total_current(device)
    top_currents.append(tot_top)
    bot_currents.append(tot_bot)
    print(f"Total top: {tot_top:.2E} and total bottom: {tot_bot:.2E}")
    # bias_volt += 0.1
    # v += 1
# plt.plot(volt_sweep, top_currents)
# plt.plot(volt_sweep, bot_currents)
plt.semilogy(volt_sweep, np.abs(np.array(top_currents) - np.array(bot_currents)))
plt.show()

write_devices(file="boring diode.dat", type="tecplot")


x=get_node_model_values(device=device, region=region, name="x")
ymax = 10
ymin = 10
fields = ("Electrons", "Holes", "Donors", "Acceptors")
for i in fields:
   node_m_vals=get_node_model_values(device=device, region=region, name=i)
   if (max(node_m_vals) > ymax):
       ymax = max(node_m_vals)
   plt.semilogy(x, node_m_vals)
plt.xlabel('x (cm)')
plt.ylabel('Density (#/cm^3)')
plt.legend(fields)
ymax *= 10
plt.axis([min(x), max(x), ymin, ymax])
plt.savefig("diode_1d_density.png")

plt.clf()
edge_average_model(device=device, region=region, node_model="x", edge_model="xmid")
xmid=get_edge_model_values(device=device, region=region, name="xmid")
efields = ("ElectronCurrent", "HoleCurrent", )
node_m_vals=get_edge_model_values(device=device, region=region, name="ElectronCurrent")
ymin=min(node_m_vals)
ymax=max(node_m_vals)
for i in efields:
 node_m_vals=get_edge_model_values(device=device, region=region, name=i)
 if min(node_m_vals) < ymin:
   ymin = min(node_m_vals)
 elif max(node_m_vals) > ymax:
   ymax = max(node_m_vals)
plt.plot(xmid, node_m_vals)
plt.xlabel('x (cm)')
plt.ylabel('J (A/cm^2)')
plt.legend(efields)
# plt.axis([min(x), max(x), 0.5*ymin, 2*ymax])
plt.savefig("diode_1d_current.png")
plt.show()
print (ymin)
print (ymax)
