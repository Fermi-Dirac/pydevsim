"""
This module holds classes and functions for dealing with Regions.
These are segments of a device which consist of one material and various parameters
These are associated with a mesh object which contains details on the simulation.

Each region is also bounded by Interfaces
"""
from pydevsim import setup_logger
from ds import set_parameter
from .material import Silicon


logger = setup_logger(__name__)
class Region(object):
    """
    Base class region
    Regions are sections within a mesh and device which have the same set of variables, or material parameters
    A 'region' generally denotes a portion of the device. Like the p-doped region, or the barrier region.
    In this way, 'regions' can be thought of as local materials with different material properties from the base Material

    We keep a 'next_id' to ensure that all region names are always unique
    """
    _next_id = 0

    def __init__(self, name='default', material=None):
        self.id = self._next_id
        self._next_id += 1
        self.name = f"{name}_id_{self.id:d}"
        if material is None:
            logger.warning(f"No material specified for region {name}, using Silicon")
            material = Silicon()
        self.material = material


class Region1D(Region):
    """
    Spceific 1D region
    Just like the base Region, but confined to 1D
    As such it only has a 'thickness', and two interfaces.

    Has a start, and end. Only has two interfaces left and right.
    """
    def __init__(self, thickness=1, left_interface=None, right_interface=None, mesh_data=None, **kwargs):
        """

        :param thickness: Thickness of this region in mesh.scale units
        :param left_interface: Interface object immediately preceeding this region
        :param right_interface: Interface object immediately following this region
        :param mesh_data: Either a list of dictionary of kwargs to be unpacked for add_line, or a list of tuples to be unpacked
                          Format matches that of Mesh1D object's add_line
        :param kwargs:
        """
        super().__init__(**kwargs)
        self.thickness = thickness
        self.left_interface = left_interface
        self.right_interface = right_interface
        if mesh_data is None:
            mesh_data = [{'position':0, 'spacing': thickness/100, 'tag': self.name + '_top'},
                         {'position': thickness, 'spacing' : thickness/100, 'tag': self.name + '_bottom'}]
        else:
            try:
                mesh_data['tag']
            except TypeError:
                logger.error("Tried to pass tuple list as mesh_data instead of dictionary. Recasting")
                mesh_data = [{'position': minfo[0], 'spacing':minfo[1], 'tag':minfo[2]} for minfo in mesh_data]
        self.mesh_data = mesh_data



