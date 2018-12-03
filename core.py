import ds

def get_ds_status(short=True):
    """
    Prints the status of the current devsim setup and all variables, solutions, node models and edge models
    :return:
    """
    for device in  ds.get_device_list():
        print("Device: " + device)
        for region in ds.get_region_list(device=device):
            print("\tRegion :" + region)
            params = ds.get_parameter_list(device=device, region=region)
            for param in params:
                val = ds.get_parameter(device=device, region=region, name=param)
                print(f"\t\t{param} = {val}")
            for node_model in ds.get_node_model_list(device=device, region=region):
                nmvals = ds.get_node_model_values(device=device, region=region, name=node_model)
                nmstr = ','.join([f'{val:.3g}' for val in nmvals])
                print(f"\t\tNode Model '{node_model}' = {nmstr!s}")
            e_models = ds.get_edge_model_list(device=device, region=region)
            for edge_model in e_models:
                emvals = ds.get_edge_model_values(device=device, region=region, name=edge_model)
                emstr = ','.join([f'{val:.3g}' for val in emvals])
                print(f"\t\tEdge Model '{edge_model}' = {emstr}")
        for interface in ds.get_interface_list(device=device):
            print("\tInterface: " + interface)
            for imodel in ds.get_interface_model_list(device=device, interface=interface):
                intstr = ','.join([f'{val:.3g}' for val in ds.get_interface_model_values(device=device, interface=interface, name=imodel)])
                print(f"\t\tInterface Model: '{imodel}' = {intstr}")
            for ieq in ds.get_interface_equation_list(device=device, interface=interface):
                # comm = ds.get_interface_equation_command(device=device, interface=interface, name=ieq)
                print(f"\t\tInterface Equation: {ieq}")
        contacts = ds.get_contact_list(device=device)
        for contact in contacts:
            print("\tContact : " + contact)
            c_eqs = ds.get_contact_equation_list(device=device, contact=contact)
            for ceq in c_eqs:
                print("\t\tContact Equation : " + ceq)


def get_ds_equations():
    devices = ds.get_device_list()
    for device in devices:
        print("Device: " + device)
        regions = ds.get_region_list(device=device)
        for region in regions:
            print("\tRegion :" + region)
            eqs = ds.get_equation_list(device=device, region=region)
            print(f"\t\t{eqs!s}")
            # params = ds.get_parameter_list(device=device, region=region)
            for eq in eqs:
                val = ds.get_equation_command(device=device, region=region, name=eq)
                print(f"\t\t{eq} = {val}")
        contacts = ds.get_contact_list(device=device)
        for contact in contacts:
            print("\tContact : " + contact)
            c_eqs = ds.get_contact_equation_list(device=device, contact=contact)
            for ceq in c_eqs:
                print("\t\tContact Equation : " + ceq)
