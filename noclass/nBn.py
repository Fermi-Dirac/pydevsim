# Gross attempt at making an nBn and using no code re-use methodologies at all. Because Fuck it.


# from ds import *
import ds
from python_packages.simple_physics import CreateSiliconDriftDiffusion, CreateSiliconSiliconInterface
from diode_common import InitialSolution
from python_packages.model_create import CreateSolution, CreateContactNodeModel
from matplotlib import pyplot as plt
# from . import diode_common
# import diode_common
q      = 1.6e-19 # coul
k      = 1.3806503e-23 # J/K
eps_0  = 8.85e-14 # F/cm^2
contactcharge_node="contactcharge_node"
contactcharge_edge="contactcharge_edge"
ece_name="ElectronContinuityEquation"
hce_name="HoleContinuityEquation"
# 1E-10 is machine epsilon problems. If n+p dominates, or ni dominates it reduces
celec_model = "(1e-10 + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
chole_model = "(1e-10 + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"

def set_dd_parameters(device, region, permittivity=11.1, n_i=1E10, T=300, mu_n=400, mu_p=200, taun=1E-5, taup=1E-5):
    names = ['Permittivity', 'ElectronCharge', 'n_i', 'T', 'kT', 'V_t', 'mu_n', 'mu_p', 'n1', 'p1', 'taun', 'taup']
    values = [permittivity*eps_0, q, n_i, T, k*T, k*T/q, mu_n, mu_p, n_i, n_i, taun, taup]
    for name, value in zip(names, values):
        ds.set_parameter(device=device, region=region, name=name, value=value)

def CreateMesh(device, region):
  '''
    Meshing
  '''
  ds.create_1d_mesh(mesh="dio")
  ds.add_1d_mesh_line(mesh="dio", pos=0, ps=1e-7, tag="top")
  ds.add_1d_mesh_line(mesh="dio", pos=0.5e-5, ps=1e-9, tag="mid")
  ds.add_1d_mesh_line(mesh="dio", pos=1e-5, ps=1e-7, tag="bot")
  ds.add_1d_contact  (mesh="dio", name="top", tag="top", material="metal")
  ds.add_1d_contact  (mesh="dio", name="bot", tag="bot", material="metal")
  ds.add_1d_region   (mesh="dio", material="Si", region=region, tag1="top", tag2="bot")
  ds.finalize_mesh(mesh="dio")
  ds.create_device(mesh="dio", device=device)

def SetSiliconParameters(device, region, T=300, esp_si=11.1, n_i=1E10, mu_n=400, mu_p=200, tau_n=1E-5, tau_p=1E-5):
  '''
    Sets physical parameters assuming constants
  '''
  #### TODO: make T a free parameter and T dependent parameters as models
  ds.set_parameter(device=device, region=region, name="Permittivity", value=esp_si * eps_0)
  ds.set_parameter(device=device, region=region, name="ElectronCharge", value=q)
  ds.set_parameter(device=device, region=region, name="n_i", value=n_i)
  ds.set_parameter(device=device, region=region, name="T", value=T)
  ds.set_parameter(device=device, region=region, name="kT", value=k * T)
  ds.set_parameter(device=device, region=region, name="V_t", value=k*T/q)
  ds.set_parameter(device=device, region=region, name="mu_n", value=mu_n)
  ds.set_parameter(device=device, region=region, name="mu_p", value=mu_p)
  #default SRH parameters
  ds.set_parameter(device=device, region=region, name="n1", value=n_i)
  ds.set_parameter(device=device, region=region, name="p1", value=n_i)
  ds.set_parameter(device=device, region=region, name="taun", value=tau_n)
  ds.set_parameter(device=device, region=region, name="taup", value=tau_p)

def SetParameters(device, region):
    return SetSiliconParameters(device, region)

def CreateNodeModel(device, region, model, expression):
  '''
    Creates a node model
  '''
  result=ds.node_model(device=device, region=region, name=model, equation=expression)
  if True:
    print(("NODEMODEL {d} {r} {m} \"{re}\"".format(d=device, r=region, m=model, re=result)))

#
def DriftDiffusionInitialSolution(device, region, circuit_contacts=None):
    ####
    #### drift diffusion solution variables
    ###
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")

    ####
    #### create initial guess from dc only solution
    ####
    ds.set_node_values(device=device, region=region, name="Electrons", init_from="IntrinsicElectrons")
    ds.set_node_values(device=device, region=region, name="Holes", init_from="IntrinsicHoles")

    ###
    ### Set up equations
    ###
    CreateSiliconDriftDiffusion(device, region)
    for i in ds.get_contact_list(device=device):
        if circuit_contacts and i in circuit_contacts:
            CreateSiliconDriftDiffusionAtContact(device, region, i, True)
        else:
            CreateSiliconDriftDiffusionAtContact(device, region, i)

#
def CreateSiliconDriftDiffusionAtContact(device, region, contact, is_circuit=False):
    '''
      Restrict electrons and holes to their equilibrium values
      Integrates current into circuit
    '''
    contact_electrons_model = "Electrons - ifelse(NetDoping > 0, {0}, n_i^2/{1})".format(celec_model, chole_model)
    contact_holes_model = "Holes - ifelse(NetDoping < 0, +{1}, +n_i^2/{0})".format(celec_model, chole_model)
    contact_electrons_name = "{0}nodeelectrons".format(contact)
    contact_holes_name = "{0}nodeholes".format(contact)

    CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
    # TODO: The simplification of the ifelse statement is time consuming
    #  CreateContactNodeModelDerivative(device, contact, contact_electrons_name, contact_electrons_model, "Electrons")
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_electrons_name, "Electrons"), "1")

    CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_holes_name, "Holes"), "1")

    #TODO: keyword args
    if is_circuit:
        ds.contact_equation(device=device, contact=contact, name="ElectronContinuityEquation", variable_name="Electrons",
                          node_model=contact_electrons_name,
                          edge_current_model="ElectronCurrent", circuit_node=GetContactBiasName(contact))

        ds.contact_equation(device=device, contact=contact, name="HoleContinuityEquation", variable_name="Holes",
                          node_model=contact_holes_name,
                          edge_current_model="HoleCurrent", circuit_node=GetContactBiasName(contact))

    else:
        ds.contact_equation(device=device, contact=contact, name="ElectronContinuityEquation", variable_name="Electrons",
                          node_model=contact_electrons_name,
                          edge_current_model="ElectronCurrent")

        ds.contact_equation(device=device, contact=contact, name="HoleContinuityEquation", variable_name="Holes",
                          node_model=contact_holes_name,
                          edge_current_model="HoleCurrent")
#####
# dio1
#
# Make doping a step function
# print dat to text file for viewing in grace
# verify currents analytically
# in dio2 add recombination
#

device = "nBn"
meshname = 'nBn_mesh'
n1_region = "n1_region"
b_region = "Barrier"
n2_region = "n2_region"
regions = [n1_region, b_region, n2_region]
n1_thick = 1e-5
b_thick = 0.5e-5
n2_thick = 10e-5
spacing = 5E-7
doping = 5E15
n_eg = 0.121 # eV
n_affinity = 0.052 # less than chi of InAs
b_eg = 1.636 # eV
b_affinity = 1.569 # less than chi of InAs
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

# Add some physics

for region in regions:
    set_dd_parameters(device=device, region=region, taun=1E-5, taup=1E-5)
    ds.node_model(device=device, region=region, name="Acceptors", equation="1.0e16")
    ds.node_model(device=device, region=region, name="Donors", equation="1.0e14")
    ds.node_model(device=device, region=region, name="NetDoping", equation="Donors-Acceptors")
    ds.print_node_values(device=device, region=region, name="NetDoping")
    InitialSolution(device, region)
    DriftDiffusionInitialSolution(device=device, region=region)
    CreateSiliconDriftDiffusion(device=device, region=region)

CreateSiliconSiliconInterface(device=device, interface=f"{n1_region}_to_{b_region}")
CreateSiliconSiliconInterface(device=device, interface=f"{b_region}_to_{n2_region}")

material = ds.get_material(device=device, region=n1_region)
print(ds.get_parameter_list(device=device, region=n1_region))

# Initial DC solution
ds.solve(type="dc", absolute_error=1.0, relative_error=1e-10, maximum_iterations=30)
# for region in regions:
#     DriftDiffusionInitialSolution(device, region)

###
### Drift diffusion simulation at equilibrium
###
currs = []
volts = []
bias_contact = "top_n1"
for volt in range(-20, 20, 1):
    volt = volt/10
    volts.append(volt)
    ds.set_parameter(device=device, name=f"{bias_contact}_bias", value=float(volt))
    ds.solve(type="dc", absolute_error=1e1, relative_error=1e-10, maximum_iterations=30)
    e_current = ds.get_contact_current(device=device, contact=bias_contact, equation="ElectronContinuityEquation")
    h_current = ds.get_contact_current(device=device, contact=bias_contact, equation="HoleContinuityEquation")
    current = e_current + h_current
    currs.append(current)
print(volts, currs)
plt.plot(volts, currs)
####
#### Ramp the bias to 0.5 Volts
####

ds.write_devices(file="nBn_ugly.dat", type="devsim")
plt.show()