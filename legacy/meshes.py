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
    thickness = thickness * 1E-5
    spacing = spacing * 1E-5  # convert nm to cm
    ds.create_1d_mesh(mesh=meshname)
    ds.add_1d_mesh_line(mesh=meshname, pos=0, ps=spacing, tag=top)
    ds.add_1d_mesh_line(mesh=meshname, pos=thickness, ps=spacing, tag=bottom)
    ds.add_1d_region(mesh=meshname, tag1=top, tag2=bottom, region=region, material=material)
    ds.add_1d_contact(material='metal', mesh=meshname, name=f"{top}_contact", tag=top)
    ds.add_1d_contact(material='metal', mesh=meshname, name=f"{bottom}_contact", tag=bottom)
    ds.finalize_mesh(mesh=meshname)
    return meshname, region, top, bottom

def p_n_junction():
    pass