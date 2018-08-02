from enum import Enum


class Material(object):
    """
    This is the base class for materials and regions
    This hosts common methods and attributes for other materials
    """
    def __init__(self, thickness, name='default'):
        self.thickness = thickness
        self.name = name

    def __repr__(self):
        return self.name

class Semiconductor(Material):
    """
    Class for solid-state simple semicondcutors and alloys

    Current version is very 1D, and is based off of AFORS style
    """
    def __init__(self, thickness, dielectric_const=11.9,
                 work_function=4.05, band_gap=1.124,
                 Ncond=3E19, Nval=3E19,
                 mobility_n=1107, mobility_p=424.6,
                 Na=0, Nd=0,
                 vel_electron=1E7, vel_holes=1E7,
                 band_to_band_tunnel_rate=0,
                 auger_e=0, auger_h=0,
                 defect_list=None,**kwargs):
        super().__init__(thickness, **kwargs)
        self.dielectric_const = dielectric_const
        self.work_function = work_function
        self.band_gap = band_gap
        self.Ncond = Ncond
        self.Nval = Nval
        self.mobility_n = mobility_n
        self.mobility_p = mobility_p
        self.Na = Na
        self.Nd = Nd
        self.velocity_e = vel_electron
        self.velocity_h = vel_holes
        self.bbr = band_to_band_tunnel_rate
        self.auger_e = auger_e
        self.auger_h = auger_h
        if defect_list is None:
            defect_list = []
        self.defect_list = defect_list


class Silicon(Semiconductor):
    """
    Simple class for textbook Silicon material
    """
    def __init__(self, thickness, **kwargs):
        super().__init__(thickness, **kwargs)

class Defect(object):
    """
    Defects within the semicondcutor

    These have a spatial position, and an energetic position.
    They also have things like their cross section, character, lifetime or density and more.
    """
    def __init__(self):
        pass


"""
Begin legacy 
"""

class Metals(Enum):
    generic = 'metal'
    Aluminum = 2
    Copper = 3
    Titanium = 4
    Tungsten = 5
    Cobalt = 6


class OldSilicon(Material):
    """
    Silicon material
    Look at `self.parameters` for the material parameters
    """
    __material_name = 'silicon'
    __defaults = {
        "Permittivity": 11.1 * eps_0,
        "n_i": 1e10,
        # mu_n and mu_p are specific for Silicon
        "mu_n": 400,
        "mu_p": 200,
        # default SRH parameters
        "n1": 1e10,
        "p1": 1e10,
        "taun": 1e-5,
        "taup": 1e-5,
    }

    def __init__(self, **kwargs):
        """
        Creates a new Silicon material. Provide alternate values for the
        parameters like this:

        s = Silicon(T=327, taun=1e16, taup=1.44e-6)
        """
        self.parameters = self.__defaults.copy()
        self.parameters.update(kwargs)


class Semiconductors(Enum):
    Silicon = Silicon
    Si = Silicon
    # PolySilicon = 2
    # Poly = 2
    # PolySi = 2
    # GaAs = 3
    # AlGaAs = 4
    # InGaAs = 5
    # SiGe = 6
    # InP = 7
    # Germanium = 8
    # Ge = 8
    # SiC_4H = 9
    # SiC_6H = 10
    # SiC_3C = 11


class Insulators(Enum):
    Photoresist = 1
    Oxide = 2
    SiO2 = 2
    Nitride = 3


class Impurities(Enum):
    Aluminum = 1
    Al = 1

    Antimony = 2
    Sb = 2

    Arsenic = 3
    As = 3

    Beryllium = 4
    Be = 4

    Boron = 5
    B = 5

    Carbon = 6
    C = 6

    Phosphorous = 7
    Phos = 7
    P = 7
