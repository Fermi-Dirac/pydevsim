# from ds import *
from ds import create_1d_mesh, add_1d_mesh_line, add_1d_region, add_1d_contact, finalize_mesh, create_device
# from python_packages.simple_physics import *
from ds import print_node_values, solve, get_node_model_values, edge_average_model, get_edge_model_values, get_contact_current
from python_packages.simple_physics import SetSiliconParameters, CreateNodeModel, CreateSolution, \
    CreateSiliconPotentialOnly, get_contact_list, CreateSiliconPotentialOnlyContact, set_parameter, GetContactBiasName
from ds import node_solution, edge_from_node_model, get_node_model_list
from python_packages.simple_physics import PrintCurrents, write_devices, set_node_values, CreateSiliconDriftDiffusionAtContact, CreateSiliconDriftDiffusion
import ds
from python_packages.simple_physics import ece_name, hce_name
from matplotlib import pyplot as plt
import numpy as np

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



def SetNetDoping(device, region):
    '''
      NetDoping
    '''



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

#---- Start ----
device = "MyDevice"
region = "MyRegion"

# Make a mesh for this device and region
# CreateMesh(device=device, region=region)
meshname = 'diode'
create_1d_mesh(mesh=meshname)
add_1d_mesh_line(mesh=meshname, pos=0, ps=1e-7, tag="top")
add_1d_mesh_line(mesh=meshname, pos=0.5e-5, ps=1e-9, tag="mid")
add_1d_mesh_line(mesh=meshname, pos=1e-5, ps=1e-7, tag="bot")
add_1d_contact(mesh=meshname, name="top", tag="top", material="metal")
add_1d_contact(mesh=meshname, name="bot", tag="bot", material="metal")
add_1d_region(mesh=meshname, material="Si", region=region, tag1="top", tag2="bot")
finalize_mesh(mesh=meshname)
create_device(mesh=meshname, device=device)
SetSiliconParameters(device, region, 300)
set_parameter(device=device, region=region, name="taun", value=1e-8)
set_parameter(device=device, region=region, name="taup", value=1e-8)

# SetNetDoping(device=device, region=region)
CreateNodeModel(device, region, "Acceptors", "1.0e18*step(0.5e-5-x)")
CreateNodeModel(device, region, "Donors", "1.0e18*step(x-0.5e-5)")
CreateNodeModel(device, region, "NetDoping", "Donors-Acceptors")

print_node_values(device=device, region=region, name="NetDoping")

# InitialSolution(device, region)
circuit_contacts=None
# Create Potential, Potential@n0, Potential@n1
# CreateSolution(device, region, "Potential")
name="Potential"
node_solution(name=name, device=device, region=region)
edge_from_node_model(node_model=name, device=device, region=region)

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

get_ds_status()
input("Pre solve")
# Initial DC solution
solve(type="dc", absolute_error=1.0, relative_error=1e-10, maximum_iterations=30)
get_ds_status()
input("Post Solve, Pre DD")
DriftDiffusionInitialSolution(device, region)
get_ds_status()
input("Post DD, pre DD solve")
solve(type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=30)
get_ds_status()
input("Post DD solve")
def get_total_current(device, ece_name="ElectronContinuityEquation", hce_name="HoleContinuityEquation"):
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


print(get_node_model_list(device=device))
plt.plot(xmid, node_m_vals)
plt.xlabel('x (cm)')
plt.ylabel('J (A/cm^2)')
plt.legend(efields)
# plt.axis([min(x), max(x), 0.5*ymin, 2*ymax])
plt.savefig("diode_1d_current.png")
plt.show()
print (ymin)
print (ymax)
