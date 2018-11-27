import ds
from adaptusim import setup_logger
logger = setup_logger(__name__)


def create_solution(device=None, regions=None, name=None):
    """
    Creates a variable to be solved during the simulation
    Also creates the edge models for the node solution so you have access on edges and nodes
    :param device:
    :param region:
    :param name:
    :return:
    """
    if device is None:
        device = ds.get_device_list()[0]
    if regions is None:
        regions = ds.get_region_list(device=device)
    if type(regions) is str:
        regions = [regions]
    try:
        regions[0]
    except IndexError:
        regions = [regions]
    for region in regions:
        ds.node_solution(name=name, device=device, region=region)
        ds.edge_from_node_model(node_model=name, device=device, region=region)


def create_node_and_derivatives(device, region, model, equation, deriv_list=None):
    """
    Creates derivatives with respect to all the variables of the equation named model in the device and region
    :param device:
    :param region:
    :param model:
    :param equation:
    :param variables:
    :return:
    """
    ds.node_model(device=device, region=region, name=model, equation=equation)
    if deriv_list is not None:
        for var in deriv_list:
            for node in ['n0', 'n1']:
                ds.node_model(device=device, region=region,
                              name=f"{model}:{var}@{node}",
                              equation=f'diff({equation},{var}@{node})')


def create_edge_and_derivatives(device, region, model, equation, deriv_list=None):
    """
    Creates derivatives with respect to all the variables of the equation named model in the device and region.

    These equations can ONLY contain the following things:
    Constants like 1E-2
    Parameters that are set in that device or region
    EdgeModel strings that are set in that area

    These equations CANNOT use
    NodeModel strings.

    If you want NodeModel strings, you need to call ds.edge_model_from_node
    :param device:
    :param region:
    :param model:
    :param equation:
    :param deriv_list:
    :return:
    """
    ds.edge_model(device=device, region=region, name=model, equation=equation)
    if deriv_list is not None:
        for var in deriv_list:
            for node in ['n0', 'n1']:
                ds.edge_model(device=device, region=region,
                              name=f'{model}:{var}@{node}',
                              equation=f'diff({equation}, {var}@{node})')
                # ds.edge_model(device=device, region=region, name=f"{model}:{var}",
                #               equation=f'simplify(diff({equation},{var}))')


def set_parameters(device=None, region=None, **kwargs):
    # TODO check if set_parameters works wihtout setting a device? Docs says it does...
    # if device is None:
    #     device = ds.get_device_list()[0]
    if device is None and region is None:
        for name, value in kwargs.items():
            ds.set_parameter(name=name, value=value)
    elif region is None:
        for name, value in kwargs.items():
            ds.set_parameter(device=device, name=name, value=value)
    else:
        for name, value in kwargs.items():
            ds.set_parameter(device=device, region=region, name=name, value=value)


def get_ds_status(short=True):
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
                if short:
                    nmvals = list(nmvals[:5]) + list(nmvals[-5:])
                print(f"\t\t Node Model '{node_model}' = {nmvals!s}")
            e_models = ds.get_edge_model_list(device=device, region=region)
            for edge_model in e_models:
                emvals = ds.get_edge_model_values(device=device, region=region, name=edge_model)
                if short:
                    emvals = emvals[0:5] + emvals[-5:0]
                print(f"\t\t Edge Model '{edge_model}' = {emvals!s}")
        contacts = ds.get_contact_list(device=device)
        for contact in contacts:
            print("\tContact : " + contact)
            c_eqs = ds.get_contact_equation_list(device=device, contact=contact)
            for ceq in c_eqs:
                print("\t\tContact Equation : " + ceq)
