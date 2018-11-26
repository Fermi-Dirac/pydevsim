import ds

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
