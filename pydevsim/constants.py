# Permitivity constants
eps_0  = 8.854187817e-14 # F/cm^2
# The electron charge quanta
Q      = 1.60217662e-19 # Couloumbs
# Boltzmann's constant
K      = 1.3806503e-23 # J/K
# Planck's constant
PLANCK_H = 6.626E-34 #J*s
# Default Ambient temperature
T      = 300 # K

ECE_NAME = "ElectronContinuityEquation"
HCE_NAME = "HoleContinuityEquation"
CELEC_MODEL = "(1e-10 + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
CHOLE_MODEL = "(1e-10 + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"