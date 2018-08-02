from os.path import join as joinpath
from os.path import dirname, abspath
from ds import set_parameter
import logging

DS_DATADIR = abspath(joinpath(dirname(__file__), 'data'))

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

log_dict = {'pydevsim.constants': logging.DEBUG,
            'pydevsim.device' : logging.DEBUG,
            'pydevsim.materials': logging.DEBUG,
            'pydevsim.mesh': logging.DEBUG}

class ParameterEnum(object):
    def set_parameters_for(self, device, region):
        """
        Use this function to register the parameters of this class into the
        simulation context. Internally uses ds.set_parameters()
        """
        props = [
            pname for pname in dir(self)
            if not pname.startswith('_') and not callable(getattr(self, pname))
        ]

        for propname in props:
            set_parameter(
                device=device,
                region=region,
                name=propname,
                value=getattr(self, propname)
            )


class _PhysicalConstants(ParameterEnum):
    # Vacuum permittivity or Dielectric constant (F/cm^2)
    eps_0 = 8.85e-14
    # The electron charge (Couloumbs)
    q = 1.6e-19
    # Planck's constant (J/K)
    k = 1.3806503e-23

PhysicalConstants = _PhysicalConstants()


class _AmbientConditions(ParameterEnum):
    # Ambient temperature (k)
    T = 300

AmbientConditions = _AmbientConditions()



def setup_logger(name, loglevel=None):
    """
    Sets up a basic logger for any of our modules
    :param name: name is normally __name__ from the module, but you could set up multiple loggers in a module if you want.
    :param loglevel: loglevel is a ENUM from logging. None implies to use the defaults above. If its not in the defaults, we do DEBUG
    :return: a logger object that lets you call its methods for setting logs. aka usage:
            mylogger = setup_logger(__name__)
            mylogger.debug("just checking")
            mylogger.critical("OMG burn it with fire!")
    """
    logger = logging.getLogger(name)
    if loglevel is None:
        loglevel = log_dict.get(name, logging.DEBUG)
    logger.setLevel(loglevel)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(loglevel)
    formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger
