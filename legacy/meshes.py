import ds

def simple_mesh(meshname='default_mesh', region='default_region', material='Si',
                top='top', bottom='bottom', thickness=10, spacing=0.1):
    """
    Creates a simple 1D mesh made of Si with metal contacts


    :param meshname:
    :param region:
    :param material:
    :param top:
    :param bottom:
    :param thickness: in nm, specify the thickness of the resistor
    :param spacing:
    :return:
    """
    thickness = thickness * 1E-4
    spacing = spacing * 1E-4  # convert nm to cm
    ds.create_1d_mesh(mesh=meshname)
    ds.add_1d_mesh_line(mesh=meshname, pos=0, ps=spacing, tag=top)
    ds.add_1d_mesh_line(mesh=meshname, pos=thickness, ps=spacing, tag=bottom)
    ds.add_1d_region(mesh=meshname, tag1=top, tag2=bottom, region=region, material=material)
    ds.add_1d_contact(material='metal', mesh=meshname, name=f"{top}_contact", tag=top)
    ds.add_1d_contact(material='metal', mesh=meshname, name=f"{bottom}_contact", tag=bottom)
    ds.finalize_mesh(mesh=meshname)
    return meshname, region, top, bottom

def p_n_junction(meshname='default_mesh', region1='p_region', region2='n_region', material1='CdS', material2='CdTe', thickness1=0.2, thickness2=2, spacing=0.05, top='top', interface='interface', bottom='bottom'):
    thickness1 = thickness1 * 1E-4
    thickness2 = thickness2 * 1E-4
    spacing = spacing * 1E-4  # convert nm to cm
    ds.create_1d_mesh(mesh=meshname)
    ds.add_1d_mesh_line(mesh=meshname, pos=0, ps=spacing/10, tag=top)
    ds.add_1d_mesh_line(mesh=meshname, pos=thickness1, ps=spacing/50, tag=interface)
    ds.add_1d_region(mesh=meshname, tag1=top, tag2=interface, region=region1, material='Si')
    ds.add_1d_mesh_line(mesh=meshname, pos=thickness1+thickness2, ps=spacing, tag=bottom)
    ds.add_1d_region(mesh=meshname, tag1=interface, tag2=bottom, region=region2, material='InAs')
    ds.add_1d_contact(material='metal', mesh=meshname, name=f"{top}_contact", tag=top)
    ds.add_1d_contact(material='metal', mesh=meshname, name=f"{bottom}_contact", tag=bottom)
    ds.add_1d_interface(mesh=meshname, tag=interface, name=interface)
    ds.finalize_mesh(mesh=meshname)
    return meshname