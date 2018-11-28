from ds import set_parameter
from os.path import abspath, join, dirname
from adaptusim import setup_logger
logger = setup_logger(__name__)
# Fundamental Constants
# Permitivity of free space
eps_0  = 8.854187817e-14 # F/cm^2
# The electron charge quanta
q      = 1.60217662e-19 # Couloumbs
# Boltzmann's constant
kb      = 1.3806503e-23 # J/K
# Planck's constant
planck_h = 6.626E-34 #J*s
# Default Ambient temperature
T      = 300 # K
# Thermal energy
kT = kb*T
# Thermal voltage
Vt = kb*T/q
# Set them globally
for name in ['eps_0', 'q', 'kb', 'planck_h', 'T', 'kT', 'Vt']:
    value = eval(name)
    logger.info(f"Now setting global variable {name} to {value:.3g}")
    set_parameter(name=name, value=value)

DS_DATADIR = abspath(join(dirname(__file__), 'data'))

ECE_NAME = "ElectronContinuityEquation"
HCE_NAME = "HoleContinuityEquation"
CELEC_MODEL = "(1e-10 + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
CHOLE_MODEL = "(1e-10 + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"

