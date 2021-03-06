"""
Mesh objects and routines
"""

import uuid
import enum
from ds import add_1d_region, add_1d_mesh_line, add_1d_interface, add_1d_contact, create_1d_mesh, finalize_mesh
from .region import Region, Region1D
from pydevsim import setup_logger
logger = setup_logger(__name__)

# class Region(object):
#     """A Region within a Semiconductor"""
#     def __init__(self, name, material):
#         super(Region, self).__init__()
#         self.name = name
#         self.material = material
#
#     def __str__(self):
#         return self.name


class Mesh(object):
    """
    The Mesh is the base of the simulation. As of now, only 1d meshes are
    implemented

    name : give a name to this Mesh. If you don't give one, a uuid will be
           used instead

    """
    _next_id = 0

    def __init__(self, name='Default_mesh', scale=1E-6):
        logger.debug(f"Creating mesh {name}")
        self.name = f"{name}_{self._next_id:d}"
        self._next_id += 1
        self.contacts = []
        self.regions = []
        self.scale = scale

    def __str__(self):
        return "Mesh id: {}\n".format(self.name) \
               + "Contacts: {}\n".format(','.join(self.contacts)) \
               + "Regions: {}\n".format(','.join(self.regions))
        # print("Mesh id: {}\n".format(self.name))
        # print("Contacts: {}\n".format(','.join(self.contacts)))
        # print("Regions: {}\n".format(','.join(self.regions)))

    def add_line(self, position, spacing, tag):
        raise NotImplementedError


    def add_contact(self, name, tag, material):
        """
        Add a contact to this mesh

        name: a name for this contact
        tag: a tag for ...?
        material: The material for the contact ...???
        """
        add_1d_contact(mesh=self.name, name=name, tag=tag, material=str(material))
        self.contacts.append(name)

    def add_region(self, name, material, tag1, tag2):
        """
        """
        raise NotImplementedError

    def finalize(self):
        """
        Last step before simulation. If you use this with a Device, it will
        call this for you.
        """
        logger.info(f"Finalizing mesh for {self.name}")
        finalize_mesh(mesh=self.name)

class Mesh1D(Mesh):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.mesh = create_1d_mesh(mesh=self.name)

    def add_line(self, position=0, spacing=0.1, tag=f"default_line"):
        """
        Add a line to the mesh

        position: Position
        spacing: Positive spacing at this location
        tag: a tag for matching contacts

        The spacing option is *positive spacing* to the next line that you have
        added.

        Additional lines are spaced from the first line using ps to the next
        specified line. The spacing of the lines are gradually spaced out to
        meet the spacing line of the next added line.

        So if you have thre lines::

            x=0.0 spacing=0.1
            x=1.0 spacing=1.0
            x=10.0 spacing=1.0

        The positions of the nodes should be about: 10 nodes from `0.0 um` to
        `0.1 um` and 10 nodes from `1.0 um` to `10.0 um` and no nodes from 10
        up.

        The purpose of this is to have tight mesh spacing in areas where it is
        important. If you had tight spacing everywhere, it could make the
        simulation much slower.
        """
        add_1d_mesh_line(
            mesh=self.name,
            pos=position * self.scale,
            ps=spacing * self.scale,
            tag=tag
        )

    def add_region(self, name, material, tag1, tag2):
        """
        name: The name of the region
        material: Any material from the enums Metals, Semiconductors

        Example:
        mesh.add_region(
            name='Cell Substrate',
            material=materials.Semiconductors.Silicon,
            tag1='top', tag2='bottom'
        )

        Example with custom parameters for Silicon:

        mesh.add_region(
            name='Cell Substrate',
            material=materials.Silicon(T=327, taun=1e16, taup=1.44e-6),
            tag1='top', tag2='bottom'
        )
        """
        # Use the material class name as material name parameter to add_1d_region?
        _mat = isinstance(material, enum.Enum) and material.__name__ or material.__class__.__name__
        add_1d_region(mesh=self.name, material=str(_mat), region=name, tag1=tag1, tag2=tag2)
        # Maybe add name, material as tuple?
        self.regions.append(Region1D(name, material))
