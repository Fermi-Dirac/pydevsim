# Gross attempt at making an nBn and using no code re-use methodologies at all. Because Fuck it.


# from ds import *
import ds
from matplotlib import pyplot as plt
q      = 1.6e-19 # coul
k      = 1.3806503e-23 # J/K
eps_0  = 8.85e-14 # F/cm^2

def set_dd_parameters(device, region, permittivity=11.1, n_i=1E10, T=300, mu_n=400, mu_p=200, taun=1E-5, taup=1E-5):
    names = ['permittivity', 'q', 'n_i', 'T', 'k', 'kT', 'V_t', 'mobility_n', 'mobility_p', 'n1', 'p1', 'taun', 'taup']
    values = [permittivity*eps_0, q, n_i, T, k, k*T, k*T/q, mu_n, mu_p, n_i, n_i, taun, taup]
    for name, value in zip(names, values):
        ds.set_parameter(device=device, region=region, name=name, value=value)

    # Setup the solutions for the potential, electron and hole densities
    pot = 'potential'
    electron_density = 'electron_density'
    hole_density = 'hole_density'
    for var in [pot, electron_density, hole_density]:
        create_solution(device=device, region=region, name=var)

    # Now for some poisson's equation
    # Create some nodemodels
    n_ie = 'n_ie', f"n_i*exp(q*{pot}/k*T)"  # Intrinsic electron density
    n_ih = 'n_ih', f'n_i^2/{n_ie[0]}'  # Intrinsic hole density
    net_n_i = 'net_n_i', f'kahan4(-{n_ie[0]}, {n_ih[0]}, p_doping, -n_doping)'  # Net intrinsic charge
    net_n_i_charge = 'net_n_i_charge', f'q*{net_n_i[0]}'  #PotentialIntrinsicCharge
    for name, eq in [n_ie, n_ih, net_n_i, net_n_i_charge]:
        ds.node_model(device=device, region=region, name=name, equation=eq)
        create_derivatives(device, region,name, eq, pot)

    E_field = 'E_field', f'({pot}@n0-{pot}@n1)*EdgeInverseLength'
    D_field = 'D_field', 'E_field * permittivity' #PotentialEdgeFlux ?? wtf

    # Initialize the electron and hole densities
    for carrier, init in zip([electron_density, hole_density], [n_ie, n_ih]):
        ds.set_node_values(device=device, region=region, name=carrier, init_from=init[0])

    # setup edge nodes
    for edge_name, eq in [E_field, D_field]:
        ds.edge_model(device=device, region=region, name=edge_name, equation=eq)
        create_derivatives(device, region, edge_name, eq, pot)

    # Create PE
    poisson_RHS= 'PoissonRHS', f'({hole_density} - {electron_density} + p_doping - n_doping)'  #*q/(permittivity)'
    # AKA pne
    ds.node_model(device=device, region=region, name=poisson_RHS[0], equation=poisson_RHS[1])
    create_derivatives(device, region, poisson_RHS[0], poisson_RHS[1], hole_density, electron_density)
    ds.equation(device=device, region=region, name="PoissonEquation", variable_name=pot,
                node_model=poisson_RHS[0], edge_model=D_field[0])

    # Use stupid bernouli formalism
    # Check if exists?
    ds.edge_from_node_model(device=device, region=region,node_model=pot)
    beta_vdiff = 'beta_vdiff', f'({pot}@n0 - {pot}@n1)*q/kT'
    ds.edge_model(device=device, region=region, name=beta_vdiff[0], equation=beta_vdiff[1])
    bernouli = 'bern', f'B({beta_vdiff[0]})'
    ds.edge_model(device=device, region=region, name = bernouli[0], equation=bernouli[1])

    # Create continuity equations
    E_qf_n =('quasi_fermi_n', f'q*({pot}-electron_affinity) + k*T*log({electron_density}/N_cond)')
    E_qf_p = ('quasi_fermi_p',f'q*({pot}-electron_affinity) - k*T*log({hole_density}/N_val) - band_gap')
    J_e = 'e_current', f'q*mobility_n*{electron_density}*EdgeInverseLength*({E_qf_n[0]}@n0-{E_qf_n[0]}@n1)'
    J_h = 'h_current', f'q*mobility_p*{hole_density}*EdgeInverseLength*({E_qf_p[0]}@n0-{E_qf_p[0]}@n1)'

    for J in [E_qf_n, E_qf_p, J_e, J_h]:
        ds.edge_model(device=device, region=region, name=J[0], equation=J[1])
        ds.node_model(device=device, region=region, name=J[0], equation=J[1])
        for node_model in [pot, electron_density, hole_density]:
            # Check if exists?!
            ds.edge_from_node_model(device=device, region=region, node_model=node_model)
            create_derivatives(device, region, J[0], J[1], node_model)

    for node_model in [n_ie, n_ih, net_n_i, net_n_i_charge]:
        ds.print_node_values(device=device, region=region, name=node_model[0])

    for edge_model in [E_qf_n, E_qf_p, J_e, J_h]:
        ds.print_edge_values(device=device, region=region, name=edge_model[0])

def create_derivatives(device, region, model, equation, *variables):
    """
    Creates derivatives with respect to all the variables of the equation named model in the device and region
    :param device:
    :param region:
    :param model:
    :param equation:
    :param variables:
    :return:
    """
    for var in variables:
        ds.node_model(device=device, region=region, name=f"{model}:{var}",
                      equation=f'simplify(diff({equation},{var}))')

def create_solution(device, region, name):
    ds.node_solution(name=name, device=device, region=region)
    ds.edge_from_node_model(node_model=name, device=device, region=region)

# Setup internal DD using old codes
# Setup new hetero interface using old codes somehow?
# or setup internal DD using E_qfn

device = "cdte_pn"
meshname = f"{device}_mesh"
CdS_region = "CdS_region"
CdTe_region = "CdTe"
regions = [CdS_region, CdTe_region]
CdS_thickness = 1E-5 * 1 # um
CdTe_thickness = 1E-5 * 10 # um
spacing = 1E-5 * 0.1 # um
p_doping = 5E15
CdTe_params = {'permittivity' : 9.4, 'electron_affinity':4.4, 'band_gap':1.5, 'N_cond':8E17, 'N_val':1.8E19,
               'mobility_n': 320, 'mobility_p':40, 'p_doping': p_doping, 'n_doping': 0}
CdS_params = {'permittivity' : 10, 'electron_affinity':4.0, 'band_gap':2.4, 'N_cond':2.2E18, 'N_val':1.8E19,
               'mobility_n': 25, 'mobility_p':100, 'p_doping': 0, 'n_doping': 1.1E18}

# Devsim requires you create the entire mesh first include interfaces before specifying regions
# Create mesh

ds.create_1d_mesh(mesh=meshname)
# n1
top_tag= f"top_{CdS_region}"
interface_tag = f"bot_{CdS_region}"
bottom_tag = f"bot_{CdTe_region}"
ds.add_1d_mesh_line(mesh=meshname, pos=0, ps=spacing, tag=top_tag)
# add_1d_mesh_line(mesh=meshname, pos=n1_thick/2, ps=spacing, tag="mid_n1")
ds.add_1d_mesh_line(mesh=meshname, pos=CdS_thickness, ps=spacing, tag=interface_tag)
ds.add_1d_region(mesh=meshname, material="CdS", region=CdS_region, tag1=top_tag, tag2=interface_tag)
ds.add_1d_contact(mesh=meshname, name=top_tag, tag=top_tag, material="metal")
# b
# add_1d_mesh_line(mesh=meshname, pos=n1_thick, ps=spacing, tag="top_b")
ds.add_1d_mesh_line(mesh=meshname, pos=CdS_thickness+CdTe_thickness, ps=spacing, tag=bottom_tag)
ds.add_1d_region(mesh=meshname, material="CdTe", region=CdTe_region, tag1=interface_tag, tag2=bottom_tag)
ds.add_1d_interface(mesh=meshname, tag=interface_tag, name=f"{CdS_region}_to_{CdTe_region}_interface")
ds.add_1d_contact(material='metal', mesh=meshname, tag=bottom_tag, name=bottom_tag)
ds.finalize_mesh(mesh=meshname)

ds.create_device(mesh=meshname, device=device)

# Add some physics

for region, params in zip(regions, [CdS_params, CdTe_params]):
    for name, value in params.items():
        ds.set_parameter(device=device, region=region, name=name, value=value)
    set_dd_parameters(device=device, region=region)


# Initial DC solution
ds.solve(type="dc", absolute_error=1.0e10, relative_error=1e10, maximum_iterations=150)
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