import logging

log_dict = {'pydevsim.constants': logging.DEBUG,
            'pydevsim.device' : logging.DEBUG,
            'pydevsim.materials': logging.DEBUG,
            'pydevsim.mesh': logging.DEBUG}

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
    if len(logger.handlers) == 0:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(loglevel)
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    return logger


# Flatten our function calls from core module
from .core import *
from .constants import *
from . import plot
from .plot import *
from . import core

# TODO Legacy
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

# TODO legacy
class _PhysicalConstants(ParameterEnum):
    # Vacuum permittivity or Dielectric constant (F/cm^2)
    eps_0 = 8.85e-14
    # The electron charge (Couloumbs)
    q = 1.6e-19
    # Planck's constant (J/K)
    k = 1.3806503e-23

PhysicalConstants = _PhysicalConstants()

# TODO legacy
class _AmbientConditions(ParameterEnum):
    # Ambient temperature (k)
    T = 300

AmbientConditions = _AmbientConditions()



