import ds
from adaptusim.legacy import create_node_and_derivatives, create_edge_and_derivatives, create_solution, set_parameters, get_ds_status
from adaptusim.legacy.meshes import simple_mesh
# List of names of things
# Fundamental
q = 'q' # fundamental unit of charge
k = 'k'  # Boltzmann's constant
eps_0  = 'eps_0' # The permittivity of free space in F/cm^2
Vt = 'Vt' # Thermal voltage = q/kt
# Solutions
potential='potential' # Electrostatic columnb potential found in Poisson's equation
delta_pot = 'vdiff' # Change in potental from one node to another normalized by thermal voltage Vt
qfn = 'quasi_fermi_n' # Quasi fermi level for electrons
qfp = 'quasi_fermi_p' # Quasi fermi level for holes
hole_density='hole_density' # Density of holes in charges / cm2
h_current_name = 'hole_current'
electron_density='electron_density'  # Density of electrons in charges / cm2
e_current_name = 'electron_current'  # Current density of electrons in charge / cm2 / s or Amps/cm2
# Nodes
total_charge='total_charge'  # The total density of net charges in charges/cm2
p_doping='p_doping' # Ionized doping concentration in charges / cm2
n_doping='n_doping' # Ionized doping concentration in charges / cm2
intrinsic_charge = 'intr_chrg'  # Total charge in columbs that is net and intrinsic
intrinsic_electrons = 'intrinsic_electrons'
intrinsic_holes = 'intrinsic_holes'
U_SRH = 'U_SRH'  # Recombination and Generation potential for Shockley-reed-hall effects
e_gen = 'e_gen'  # Generation rate of electrons
h_gen = 'h_gen'  # Generation rate for holes
# Edges
e_field = 'e_field'  # The electric field approximated by delta potential / distance
d_field = 'd_field'  # The displacement field which is the electric field damped by the permittivity in a material
# Parameters
permittivity = 'permittivity'  # the dielectric coefficient of a material or free space

# End of names

def setup_physical_constants(device=None, region=None, q=1.6E-19, k=1.38065E-23, eps_0=8.85E-14, **kwargs):
    return set_parameters(**locals())


def setup_poisson_parameters(device, region, T=300, permittivity=11.1, n_i=1E10,
                             mobility_n=400, mobility_p=200,
                             p_doping=1E12, n_doping=1E16,
                             tau_n=1E-5, tau_p=1E-5):
    ds.node_model(device=device, region=region, name=Vt, equation='k*T/q')
    return set_parameters(**locals())


def setup_poisson(device, region, variable_update='log_damp', **kwargs):
    """
    Sets up poisson equation solution
    Requires:
     setup physical constants and setup poisson parameters
    :param kwargs:
    :return:
    """
    # These are variables that are to be solved during the simulation
    # for sol_variable in [hole_density, electron_density, potential, qfn, qfp]:
    for sol_variable in [potential, qfn, qfp]:
        create_solution(device, region, sol_variable)
    hole_density_eq = f'n_i*exp( ({qfp}-{potential})*k*T/{q})'
    electron_density_eq = f'n_i*exp(({potential}-{qfn})*k*T/{q})'
    # ds.equation(device, region, name='HoleQFL', variable_name=qfn,
    #             node_model=hole_density, edge_model=)
    total_charge_eq = f'-{q}*kahan4({hole_density}, -{electron_density}, {p_doping}, -{n_doping})'

    for name, eq in zip([electron_density, hole_density, total_charge],
                        [electron_density_eq, hole_density_eq, total_charge_eq]):
        create_node_and_derivatives(device, region, name, eq, [potential])

    # # These are variables that exist per each node and may be spatially dependant
    # intrinsic_charge_eq = f'-{q}*kahan4({intrinsic_holes}, -{intrinsic_electrons}, {p_doping}, -{n_doping})'
    # create_node_and_derivatives(device, region, total_charge, total_charge_eq, [hole_density, electron_density, potential])
    # # intrinsic_electrons_eq = f'n_i*exp(({potential}-{qfn})/(k*T/q))'  # TODO this is not intrinsic, thsi is any!
    # # intrinsic_holes_eq = f'n_i*exp(({qfp}-{potential})/(k*T/q))'
    # # Homojunction
    # intrinsic_holes_eq = f'{10E10}^2/{intrinsic_electrons}'# f'{intrinsic_charge}^2/{intrinsic_electrons}'
    # intrinsic_electrons_eq = f"{Vt}*exp({potential}*{Vt})"# f"{intrinsic_charge}*exp({potential}*q/(k*T))"
    #
    # for name, eq in zip([intrinsic_electrons, intrinsic_holes, intrinsic_charge],
    #                     [intrinsic_electrons_eq, intrinsic_holes_eq, intrinsic_charge_eq]):
    #     create_node_and_derivatives(device, region, name, eq, [potential])

    # for name, eq in zip([qfn, qfp], [])
    # These are the variables for each edge, and connect nodes together
    e_field_eq = f"({potential}@n0-{potential}@n1)*EdgeInverseLength"
    d_field_eq = f'{e_field}*{permittivity}*{eps_0}'
    for name, eq in zip([e_field, d_field],[e_field_eq, d_field_eq]):
        create_edge_and_derivatives(device, region, name, eq, [potential])

    # Lastly, setup the equation to solve for 'potential'
    # ds.set_node_values(device=device, region=region, name=qfp, init_from=intrinsic_holes)
    # ds.set_node_values(device=device, region=region, name=qfn, init_from=intrinsic_electrons)
    ds.equation(device=device, region=region, name='PoissonEq', variable_name=potential,
                node_model=total_charge, edge_model=d_field, variable_update=variable_update)
    # ds.set_node_values(device=device, region=region, name=potential, init_from='0')


def setup_continuity(device, region, **kwargs):
    """

    :param device:
    :param region:
    :param kwargs:
    :return:
    """
    # These are variables that are to be solved during the simulation
    # TODO check if these already exist?
    for sol_variable in [hole_density, electron_density, potential]:
        create_solution(device, region, sol_variable)

    # These are variables that exist per each node and may be spatially dependant
    # Generation rate
    U_SRH_eq = f'{q}*({electron_density}*{hole_density} - {intrinsic_charge}^2)/(tau_p*({electron_density}+{intrinsic_electrons}) + tau_n*({hole_density}+{intrinsic_holes}))'
    for name, eq in zip([U_SRH, e_gen, h_gen], [U_SRH_eq, '-'+U_SRH_eq, U_SRH_eq]) :
        # TODO can we move SRH to its own function/method/module? yes, yes we can.
        create_node_and_derivatives(device, region, name, eq, [electron_density, hole_density])

    # Bernouli setup
    ds.node_model(device=device, region=region, name=Vt, equation='k*T/q')
    ds.node_model(device=device, region=region, name=Vt+"_inv", equation='q/(k*T)')
    for name, eq in [
        (f'{delta_pot}', f"q/(k*T)*({potential}@n0 - {potential}@n1)"),
                     (f'Bern01', f'B({delta_pot})'),
                     (f'{delta_pot}:{potential}@n0', f'q/(k*T)'),
                     (f'{delta_pot}:{potential}@n1', f'-{delta_pot}:{potential}@n0'),
                     (f'Bern01:{potential}@n0', f'dBdx({delta_pot})*{delta_pot}:{potential}@n0'),
                     (f'Bern01:{potential}@n1', f'-Bern01:{potential}@n0')]:
        ds.edge_model(device=device, region=region, name=name, equation=eq)
    # electron continuity
    e_current = f"{q}*mobility_n*EdgeInverseLength*k*T/{q}*kahan3({electron_density}@n1*Bern01, {electron_density}@n1*{delta_pot}, -{electron_density}@n0*Bern01)"
    create_edge_and_derivatives(device, region, e_current_name, e_current, [electron_density, hole_density,potential])
    ds.equation(device=device, region=region, name='ECE', variable_name=electron_density,
                time_node_model=electron_density,
                edge_model=e_current_name, node_model=e_gen)

    # hole continuity
    h_current = f"-{q}*mobility_p*EdgeInverseLength*k*T/{q}*kahan3({hole_density}@n1*Bern01, -{hole_density}@n0*vdiff, -{hole_density}@n0*Bern01)"
    create_edge_and_derivatives(device, region, h_current_name, h_current, [electron_density, hole_density, potential])
    ds.equation(device=device, region=region, name='HCE', variable_name=hole_density,
                time_node_model=hole_density,
                edge_model=h_current_name, node_model=h_gen)
    # ds.equation(device=device, region=region, name=, variable_name=electron_density,
    #             edge_model=e_current, node_model=e_gen)


    # ds.equation(device=device, region=region, name=, variable_name=hole_density,
    #             edge_model=, node_model=h_gen)
    # These are edge variables between nodes


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    device = 'testdevice'
    region = 'testregion'
    meshname = 'testmesh'
    simple_mesh(meshname=meshname, region=region)
    ds.create_device(mesh=meshname, device=device)
    setup_physical_constants(device, region)
    setup_poisson_parameters(device, region, p_doping=1e18)
    setup_poisson(device, region)
    get_ds_status()
    ds.solve(type="dc", absolute_error=10, relative_error=10, maximum_iterations=30)
    get_ds_status()
    print("Post solve")
    setup_continuity(device, region)
    ds.solve(type="dc", absolute_error=10, relative_error=1e1, maximum_iterations=30)
    for volt in [-1, 0, 1]:
        ds.set_parameter(device=device, name='top_contact_bias', value=float(volt))
        chrg = ds.get_node_model_values(device=device, region=region, name=total_charge)
        pot = ds.get_edge_model_values(device=device, region=region, name=e_field)
        for data in [chrg, pot]:
            plt.figure()
            plt.plot(data)
        plt.show()
# Ei = (Ec + Ev) / 2