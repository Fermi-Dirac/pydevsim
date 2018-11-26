from os.path import join as joinpath
from os.path import dirname, abspath
from ds import set_parameter
import ds
import logging
import numpy as np
try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None
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
V_t = kb*T/q
# Set them globally
for name in ['eps_0', 'q', 'kb', 'planck_h', 'T', 'kT', 'V_t']:
    value = eval(name)
    print(f"Now setting global variable {name} to {value:.3g}")
    set_parameter(name=name, value=value)

DS_DATADIR = abspath(joinpath(dirname(__file__), 'data'))

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

def plot_charge(device=None,regions=None, charge_names=None):
    """
    Plots the charge along the device
    :param device:
    :param regions:
    :param charge_names:
    :return:
    """
    if device is None:
        device = ds.get_device_list()[0]
    if regions is None:
        regions = ds.get_region_list(device=device)
    if charge_names is None:
        charge_names = ("Electrons", "Holes", "Donors", "Acceptors")
    plt.figure()
    for var_name in charge_names:
        total_x = np.array([])
        total_y = []
        for region in regions:
            x = np.array(ds.get_node_model_values(device=device, region=region, name="x"))
            # print(var_name, min(x), max(x), min(x)*1e4, max(x)*1e4)
            total_x = np.append(total_x, x)
            y=ds.get_node_model_values(device=device, region=region, name=var_name)
            total_y.extend(y)
        # plt.axis([min(x), max(x), ymin, ymax])
        plt.semilogy(np.array(total_x)*1e4, total_y)
    plt.xlabel('x (um)')
    plt.ylabel('Density (#/cm^3)')
    plt.legend(charge_names)
    # plt.savefig("diode_1d_density.png")
    # plt.show()

def plot_current(device=None, regions=None, current_names=None):
    if device is None:
        device = ds.get_device_list()[0]
    if regions is None:
        regions = ds.get_region_list(device=device)
    if current_names is None:
        current_names = ('Jn', 'Jp')
    plt.figure()
    for var_name in current_names:
        total_x = np.array([])
        total_y = []
        for region in regions:
            ds.edge_average_model(device=device, region=region, node_model="x", edge_model="xmid")
            x = np.array(ds.get_edge_model_values(device=device, region=region, name="xmid"))
            # print(var_name, min(x), max(x), min(x)*1e4, max(x)*1e4)
            total_x = np.append(total_x, x)
            y=ds.get_edge_model_values(device=device, region=region, name=var_name)
            total_y.extend(y)
        # plt.axis([min(x), max(x), ymin, ymax])
        plt.plot(np.array(total_x)*1e4, total_y)
    plt.xlabel('x (um)')
    plt.ylabel('Current (A/cm^2)')
    plt.legend(current_names)
    # plt.show()

def plot_potential(device=None, regions=None, potential='Potential'):
    if device is None:
        device = ds.get_device_list()[0]
    if regions is None:
        regions = ds.get_region_list(device=device)
    plt.figure()
    pots = np.array([])
    total_x = np.array([])
    for region in regions:
        x = np.array(ds.get_node_model_values(device=device, region=region, name="x"))
        # print(var_name, min(x), max(x), min(x)*1e4, max(x)*1e4)
        total_x = np.append(total_x, x)
        pots = np.append(pots, ds.get_node_model_values(device=device, region=region, name=potential))
    plt.plot(total_x*1e4, pots)
    plt.xlabel('X (um)')
    plt.ylabel('Potential (V)')


def get_ds_status():
    """
    Prints the status of the current devsim setup and all variables, solutions, node models and edge models
    :return:
    """
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