from matplotlib import pyplot as plt
import ds
import numpy as np
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
    plt.figure(figsize=(6, 2.7))
    labels = []
    for var_name in charge_names:
        total_x = np.array([])
        total_y = []
        for region in regions:
            labels.append(f"{var_name} in {region}")
            x = np.array(ds.get_node_model_values(device=device, region=region, name="x"))
            # print(var_name, min(x), max(x), min(x)*1e4, max(x)*1e4)
            total_x = np.append(total_x, x)
            y=ds.get_node_model_values(device=device, region=region, name=var_name)
            total_y.extend(y)
        # plt.axis([min(x), max(x), ymin, ymax])
            plt.semilogy(x*1e4, y)
    plt.xlabel('x (um)')
    plt.ylabel('Density (#/cm^3)')
    plt.legend(labels)
    plt.title('Charge Density')
    # plt.savefig("diode_1d_density.png")
    # plt.show()

def plot_current(device=None, regions=None, current_names=None, title='Current Density'):
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
    plt.title(title)
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
        new_pots = ds.get_node_model_values(device=device, region=region, name=potential)
        pots = np.append(pots, new_pots)
        plt.plot(x*1e4, new_pots, marker='o', linestyle='-', markersize=3)
    plt.legend(regions)
    # plt.plot(total_x*1e4, pots)
    plt.xlabel('X (um)')
    plt.ylabel('Potential (V)')
    plt.title(potential)

def plot_band_diagram(device=None, regions=None, ec_name='EC', ev_name='EV'):
    if device is None:
        device = ds.get_device_list()[0]
    if regions is None:
        regions = ds.get_region_list(device=device)
    plt.figure(figsize=(6, 2.7))
    labels = []
    for region in regions:
        x = np.array(ds.get_node_model_values(device=device, region=region, name="x"))
        for band in [ec_name, ev_name]:
            band_edge = ds.get_node_model_values(device=device, region=region, name=band)
            plt.plot(x * 1e4, band_edge, marker='o', linestyle='-', markersize=3)
            labels.append(f'{band} in {region}')
    plt.legend(labels)
    plt.xlabel('x (µm)')
    plt.ylabel('Energy (eV)')
    plt.title('Band Diagram')
