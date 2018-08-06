import uuid

import numpy as np
from pydevsim import setup_logger
from ds import edge_from_node_model, equation, node_model,get_node_model_list,get_edge_model_list,edge_model
from ds import contact_node_model, get_contact_current, contact_equation, set_node_values, get_contact_list, node_solution
from ds import create_device, set_parameter, solve, write_devices, get_parameter, add_1d_contact, add_1d_region
from pydevsim import ECE_NAME, HCE_NAME, CELEC_MODEL, CHOLE_MODEL # ece_name, hce_name, celec_model, chole_model
from .mesh import Mesh1D
from .environment import Environment
# from ds import *
# namespace pollution is frowned upon. Especially top-level packing imports will cause big namespace problems later

logger = setup_logger(__name__) #  logging.getLogger("Device")
#TODO: move this to the appropiate place
contactcharge_node="contactcharge_node"
contactcharge_edge="contactcharge_edge"
# ece_name="ElectronContinuityEquation"
# hce_name="HoleContinuityEquation"
# celec_model = "(1e-10 + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
# chole_model = "(1e-10 + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"



def ensure_edge_from_node_model_exists(device, region, nodemodel):
    """
    Checks if the edge models exists
    """
    if nodemodel not in get_node_model_list(device=device,region=region):
        raise ValueError(f"{nodemodel} must exist")

    # emlist = get_edge_model_list(device=device, region=region)
    emtest = ("{0}@n0".format(nodemodel) and "{0}@n1".format(nodemodel))
    if not emtest:
        logger.debug("INFO: Creating ${0}@n0 and ${0}@n1".format(nodemodel))
        edge_from_node_model(device=device, region=region, node_model=nodemodel)

# CamelCase is reserved from Class objects. Methods and functions should be snake_case (PEP8)
def create_electron_current(device, region, mu_n):
    """
    Sets up the electron current
    :param device:
    :param region:
    :param mu_n: mobility
    :return:
    """
    ensure_edge_from_node_model_exists(device, region, "Potential")
    ensure_edge_from_node_model_exists(device, region, "Electrons")
    ensure_edge_from_node_model_exists(device, region, "Holes")
    # Make sure the bernoulli functions exist
    if "Bern01" not in get_edge_model_list(device=device, region=region):
        create_bernoulli(device, region)
    # test for requisite models here
    #  Jn = "ElectronCharge*{0}*EdgeInverseLength*V_t*(Electrons@n1*Bern10 - Electrons@n0*Bern01)".format(mu_n)
    Jn = "ElectronCharge*{0}*EdgeInverseLength*V_t*kahan3(Electrons@n1*Bern01,  Electrons@n1*vdiff,  -Electrons@n0*Bern01)".format(mu_n)
    #  Jn = "ElectronCharge*{0}*EdgeInverseLength*V_t*((Electrons@n1-Electrons@n0)*Bern01 +  Electrons@n1*vdiff)".format(mu_n)
    edge_model(device=device, region=region, name="ElectronCurrent", equation=Jn)
    for i in ("Electrons", "Potential", "Holes"):
        create_edge_model_derivatives(device, region, "ElectronCurrent", Jn, i)

def create_hole_current(device, region, mu_p):
    """
    Hole current
    """
    ensure_edge_from_node_model_exists(device, region, "Potential")
    ensure_edge_from_node_model_exists(device, region, "Holes")
    # Make sure the bernoulli functions exist
    if "Bern01" not in get_edge_model_list(device=device, region=region):
        create_bernoulli(device, region)
    # test for requisite models here
    #  Jp ="-ElectronCharge*{0}*EdgeInverseLength*V_t*(Holes@n1*Bern01 - Holes@n0*Bern10)".format(mu_p)
    Jp ="-ElectronCharge*{0}*EdgeInverseLength*V_t*kahan3(Holes@n1*Bern01, -Holes@n0*Bern01, -Holes@n0*vdiff)".format(mu_p)
    #  Jp ="-ElectronCharge*{0}*EdgeInverseLength*V_t*((Holes@n1 - Holes@n0)*Bern01 - Holes@n0*vdiff)".format(mu_p)
    edge_model(device=device, region=region, name="HoleCurrent", equation=Jp)
    for i in ("Holes", "Potential", "Electrons"):
        create_edge_model_derivatives(device, region, "HoleCurrent", Jp, i)


def create_pe(device, region):
    pne = "-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)"
    create_node_model(device, region, "PotentialNodeCharge", pne)
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Electrons")
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Holes")

    equation(device=device,
        region=region,
        name="PotentialEquation", variable_name="Potential",
        node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
        time_node_model="", variable_update="log_damp")


def create_bernoulli(device, region):
    """
    Creates the Bernoulli function for Scharfetter Gummel
    """
    # test for requisite models here
    ensure_edge_from_node_model_exists(device, region, "Potential")
    vdiffstr = "(Potential@n0 - Potential@n1)/V_t"
    for model, expression in [("vdiff", "(Potential@n0 - Potential@n1)/V_t"),
                              ("vdiff:Potential@n0", "V_t^(-1)"),
                              ("vdiff:Potential@n1", "-vdiff:Potential@n0"),
                              ("Bern01","B(vdiff)", ),
                              ("Bern01:Potential@n0","dBdx(vdiff) * vdiff:Potential@n0"),
                              ("Bern01:Potential@n1","-Bern01:Potential@n0")]:
        result = edge_model(device=device, region=region, name=model, equation=expression)
        logger.debug(f"New edgemodel {device} {region} {model} -> '{result}'")#.format(d=device, r=region, m=model, re=result))


def create_SRH(device, region):
    """
    Creates a Shockley-Reed-Hall recombination equation.
    :param device:
    :param region:
    :return:
    """
    USRH = "(Electrons*Holes - n_i^2)/(taup*(Electrons + n1) + taun*(Holes + p1))"
    Gn = "-ElectronCharge * USRH"
    Gp = "+ElectronCharge * USRH"
    create_node_model(device, region, "USRH", USRH)
    create_node_model(device, region, "ElectronGeneration", Gn)
    create_node_model(device, region, "HoleGeneration", Gp)
    for i in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "USRH", USRH, i)
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, i)
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, i)


def create_electron_continuity_eq(device, region, mu_n):
    """
    Creates the electron continuity equation to be solved
    :param device:
    :param region:
    :param mu_n: mobility of electrons
    :return:
    """
    create_electron_current(device, region, mu_n)

    NCharge = "ElectronCharge * Electrons"
    create_node_model(device, region, "NCharge", NCharge)
    CreateNodeModelDerivative(device, region, "NCharge", NCharge, "Electrons")

    equation(
        device=device, region=region,
        name="ElectronContinuityEquation", variable_name="Electrons",
        time_node_model="NCharge",
        edge_model="ElectronCurrent", variable_update="positive", node_model="ElectronGeneration"
    )


def create_hole_continuity_eq(device, region, mu_p):
    """
    Creates the hole continuity equation
    :param device:
    :param region:
    :param mu_p: mobility of holes
    :return:
    """
    create_hole_current(device, region, mu_p)
    PCharge = "-ElectronCharge * Holes"
    create_node_model(device, region, "PCharge", PCharge)
    CreateNodeModelDerivative(device, region, "PCharge", PCharge, "Holes")

    equation(
        device=device, region=region, name="HoleContinuityEquation",
        variable_name="Holes",
        time_node_model="PCharge",
        edge_model="HoleCurrent", variable_update="positive",
        node_model="HoleGeneration"
    )


# TODO: Move this to the appropiate place
def create_node_model(device, region, model, expression):
    """
        Creates a node model
    """
    result = node_model(
        device=device,
        region=region,
        name=model,
        equation=expression
    )
    logger.debug("NODEMODEL {d} {r} {m} \"{re}\"".format(d=device, r=region, m=model, re=result))


# Retired
# def InNodeModelList(device, region, model):
#     """
#     Checks to see if this node model is available on device and region
#     """
#     return model in get_node_model_list(device=device, region=region)


# Retired
# def InEdgeModelList(device, region, model):
#     """
#         Checks to see if this edge model is available on device and region
#     """
#     return model in get_edge_model_list(device=device, region=region)


# TODO: Move this to the appropiate place
def CreateEdgeModel(device, region, model, expression):
    """
    Creates an edge model
    """
    result = edge_model(device=device, region=region, name=model, equation=expression)
    logger.debug("EDGEMODEL {d} {r} {m} \"{re}\"".format(d=device, r=region, m=model, re=result))


# TODO: Move this to the appropiate place
def CreateNodeModelDerivative(device, region, model, expression, *args):
    """
    Create a node model derivative
    """
    for v in args:
        create_node_model(
            device, region,
            "{m}:{v}".format(m=model, v=v),
            "simplify(diff({e},{v}))".format(e=expression, v=v)
        )


def CreateContactNodeModel(device, contact, model, expression):
    """
    Creates a contact node model
    """
    result = contact_node_model(device=device, contact=contact, name=model, equation=expression)
    logger.debug("CONTACTNODEMODEL {d} {c} {m} \"{re}\"".format(d=device, c=contact, m=model, re=result))


# TODO: Move this to the appropiate place
def CreateSiliconPotentialOnly(device, region):
    """
        Creates the physical models for a Silicon region
    """
    if not "Potential" in get_node_model_list(device=device, region=region):
        logger("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")
    # require NetDoping
    intrinsics = (
        ("IntrinsicElectrons", "n_i*exp(Potential/V_t)"),
        ("IntrinsicHoles", "n_i^2/IntrinsicElectrons"),
        ("IntrinsicCharge", "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"),
        ("PotentialIntrinsicCharge", "-ElectronCharge * IntrinsicCharge")
    )
    for name, eq in intrinsics:
        create_node_model(device, region, name, eq)
        CreateNodeModelDerivative(device, region, name, eq, "Potential")

    # TODO: Edge Average Model
    electrics = (
        ("ElectricField", "(Potential@n0-Potential@n1)*EdgeInverseLength"),
        ("PotentialEdgeFlux", "Permittivity * ElectricField")
    )
    for name, eq in electrics:
        CreateEdgeModel(device, region, name, eq)
        create_edge_model_derivatives(device, region, name, eq, "Potential")

    equation(
        device=device, region=region,
        name="PotentialEquation", variable_name="Potential",
        node_model="PotentialIntrinsicCharge",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp"
    )


# TODO: Move this to the appropiate place
def CreateSiliconPotentialOnlyContact(device, region, contact, is_circuit=False):
    """
    Creates the potential equation at the contact
    if is_circuit is true, than use node given by GetContactBiasName
    """
    # Means of determining contact charge
    # Same for all contacts
    if "contactcharge_node" not in get_node_model_list(device=device, region=region):
        create_node_model(device, region, "contactcharge_node", "ElectronCharge*IntrinsicCharge")
    # TODO: This is the same as D-Field
    if "contactcharge_edge" not in get_edge_model_list(device=device, region=region):
        CreateEdgeModel(device, region, "contactcharge_edge", "Permittivity*ElectricField")
        create_edge_model_derivatives(device, region, "contactcharge_edge", "Permittivity*ElectricField", "Potential")


# TODO: Move this to the appropiate place
def CreateSiliconDriftDiffusion(device, region, mu_n="mu_n", mu_p="mu_p"):
    create_pe(device, region)
    create_bernoulli(device, region)
    create_SRH(device, region)
    create_electron_continuity_eq(device, region, mu_n)
    create_hole_continuity_eq(device, region, mu_p)


# TODO: Move this to the appropiate place
def CreateSiliconDriftDiffusionAtContact(device, region, contact, is_circuit=False):
    """
        Restrict electrons and holes to their equilibrium values
        Integrates current into circuit
    """
    contact_electrons_model = "Electrons - ifelse(NetDoping > 0, {0}, n_i^2/{1})".format(CELEC_MODEL, CHOLE_MODEL)
    contact_holes_model = "Holes - ifelse(NetDoping < 0, +{1}, +n_i^2/{0})".format(CELEC_MODEL, CHOLE_MODEL)
    contact_electrons_name = "{0}nodeelectrons".format(contact)
    contact_holes_name = "{0}nodeholes".format(contact)

    CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
    # TODO: The simplification of the ifelse statement is time consuming
    #  CreateContactNodeModelDerivative(device, contact, contact_electrons_name, contact_electrons_model, "Electrons")
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_electrons_name, "Electrons"), "1")

    CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_holes_name, "Holes"), "1")

    # TODO: keyword args
    if is_circuit:
        contact_equation(
            device=device,
            contact=contact,
            name="ElectronContinuityEquation",
            variable_name="Electrons",
            node_model=contact_electrons_name,
            edge_current_model="ElectronCurrent",
            circuit_node=GetContactBiasName(contact)
        )

        contact_equation(
            device=device,
            contact=contact,
            name="HoleContinuityEquation",
            variable_name="Holes",
            node_model=contact_holes_name,
            edge_current_model="HoleCurrent",
            circuit_node=GetContactBiasName(contact)
        )
    else:
        contact_equation(
            device=device,
            contact=contact,
            name="ElectronContinuityEquation",
            variable_name="Electrons",
            node_model=contact_electrons_name,
            edge_current_model="ElectronCurrent")

        contact_equation(
            device=device,
            contact=contact,
            name="HoleContinuityEquation",
            variable_name="Holes",
            node_model=contact_holes_name,
            edge_current_model="HoleCurrent")


# TODO: Move this to the appropiate place
def create_edge_model_derivatives(device, region, model, expression, variable):
    """
    Creates edge model derivatives
    """
    for num in ['n0', 'n1']:
        edge_model(device=device,
                   region=region,
                   name=f"{model}:{variable}@{num}",
                   equation = f"simplify(diff({expression}, {variable}@{num}))")
    # CreateEdgeModel(
    #     device, region,
    #     "{m}:{v}@n0".format(m=model, v=variable),
    #     "simplify(diff({e}, {v}@n0))".format(e=expression, v=variable)
    # )
    # CreateEdgeModel(
    #     device, region,
    #     "{m}:{v}@n1".format(m=model, v=variable),
    #     "simplify(diff({e}, {v}@n1))".format(e=expression, v=variable)
    # )


# TODO: Move this to the appropiate place
def drift_diffusion_initial_solution(device, region, circuit_contacts=None):
    # drift diffusion solution variables
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")

    # create initial guess from dc only solution
    set_node_values(device=device, region=region, name="Electrons", init_from="IntrinsicElectrons")
    set_node_values(device=device, region=region, name="Holes", init_from="IntrinsicHoles")

    # Set up equations
    CreateSiliconDriftDiffusion(device, region)
    for i in get_contact_list(device=device):
        if circuit_contacts and i in circuit_contacts:
            CreateSiliconDriftDiffusionAtContact(device, region, i, True)
        else:
            CreateSiliconDriftDiffusionAtContact(device, region, i)


# TODO: Move this to the appropiate place
def CreateSolution(device, region, name):
    """
        Creates solution variables
        As well as their entries on each edge
    """
    node_solution(name=name, device=device, region=region)
    edge_from_node_model(node_model=name, device=device, region=region)


class Device(object):
    """docstring for Device"""

    __models = None
    _next_id = 0

    def __init__(self, name=f"default_device"):
        self.name =f"{name}_{self._next_id:d}"
        self._next_id += 1
        # if mesh is None:
        #     raise NotImplementedError
        # create_device(mesh=self.mesh.name, device=self.name)
        self._models = []

    def _contact_bias_name(self, contact):
        return "{0}_bias".format(contact)

    # def print_currents(self):
    #     for c in self.mesh.contacts:
    #         e_current = get_contact_current(
    #             device=device, contact=c, equation='ElectronContinuityEquation'
    #         )
    #         h_current = get_contact_current(
    #             device=device, contact=c, equation='HoleContinuityEquation'
    #         )
    #         total_current = e_current + h_current
    #         voltage = get_parameter(
    #             device=device,
    #             name=self._contact_bias_name(contact)
    #         )
    #     logger.info("{0}\t{1}\t{2}\t{3}\t{4}".format(contact, voltage, e_current, h_current, total_current))

    def create_node_model(self, region, model, expression):
        result = node_model(device=self.name, region=region, name=model, equation=expression)
        logger.debug("NODEMODEL {d} {r} {m} \"{re}\"".format(d=self.name, r=region, m=model, re=result))

    def create_solution(self, region, name):
        '''
        Creates solution variables
        As well as their entries on each edge
        '''
        node_solution(name=name, device=self.name, region=region)
        edge_from_node_model(node_model=name, device=self.name, region=region)

    def solve(self, *args, **kwargs):
        if not args and not kwargs:
            self.initial_solution()
        elif kwargs.get('type', '') == 'ramp':
            kwargs['type'] = 'dc'
            start = kwargs.pop('start')
            stop = kwargs.pop('stop')
            step = kwargs.pop('step')
            contact = kwargs.pop('contact')
            for v in range(start, stop, step):
                set_parameter(
                    device=self.name,
                    name=self._contact_bias_name(contact),
                    value=v
                )
                solve(**kwargs)
                self.print_currents()
        else:
            solve(*args, **kwargs)

        for model in self.__models:
            model.solve(*args, **kwargs)

    def initial_solution(self):
        self.setup_context()
        for region in self.mesh.regions:
            # Create Potential, Potential@n0, Potential@n1
            CreateSolution(self.name, region, "Potential")

            # Create potential only physical models
            # TODO: move to materials, relate region with material
            CreateSiliconPotentialOnly(self.name, region)

            # Set up the contacts applying a bias
            # TODO: Try to use self.contacts instead
            # it is more correct for the bias to be 0, and it looks like there is side effects

            for c in self.mesh.contacts:
                set_parameter(device=self.name, name=self._contact_bias_name(c), value=0.0)
                # TODO: move to models module
                CreateSiliconPotentialOnlyContact(self.name, region, c)

    def create_solution_variable(self, name):
        """
        Creates solution variables
        As well as their entries on each edge
        """
        for region in self.mesh.regions:
            node_solution(name=name, device=self.name, region=region)
            edge_from_node_model(node_model=name, device=self.name, region=region)

    def drift_diffusion_initial_solution(self):
        # TODO: move it to somewhere else
        # drift diffusion solution variables
        self.create_solution_variable("Electrons")
        self.create_solution_variable("Holes")

        for region in self.mesh.regions:
            # Create initial guess from DC only solution
            set_node_values(
                device=self.name,
                region=region.name,
                name="Electrons",
                init_from="IntrinsicElectrons"
            )
            set_node_values(
                device=self.name,
                region=region.name,
                name="Holes",
                init_from="IntrinsicHoles"
            )

            # Set up equations
            CreateSiliconDriftDiffusion(self.name, region.name)
            for c in self.mesh.contacts:
                CreateSiliconDriftDiffusionAtContact(self.name, region.name, c)

    # TODO: This thing should be set up by the material
    def setup_context(self):
        """
            Initialize the context for the equations. Basically sets some variables.
            Region is an instance of the mesh.Region class
        """
        # Region context
        from pydevsim import T, K, Q
        for region in self.mesh.regions:
            for n, v in region.material.parameters.items():
                set_parameter(device=self.name, region=region.name, name=n, value=v)
            # General
            set_parameter(device=self.name, region=region.name, name="ElectronCharge", value=Q)
            set_parameter(device=self.name, region=region.name, name="T", value=T)
            set_parameter(device=self.name, region=region.name, name="kT", value=K*T)
            set_parameter(device=self.name, region=region.name, name="V_t", value=K*T/Q)

    def setup_model(self, model):
        self.__models.append(model)

    def export(self, filename, format='devsim_data'):
        write_devices(file=filename, type=format)


class Device1D(Device):
    """
    Simple 1D device

    """
    def __init__(self, regions, mesh=Mesh1D(), interfaces=None, contacts=None, environment=Environment(), R_series=0, **kwargs):
        """
        Creates a 1D device

        This constructor is responsbile for:
        Creating the mesh
        * Adding the meshlines from the Region objects inside of regions
        * Adding interfaces between the Regions
        * Upading the Region objects interfaces as needed
        * Creating the contacts as specified by giving a list of two Regions
        :param regions: List of Region objects which create this 1D device. Order matters and is from top to bottom
        :param interfaces: List of interfaces for this device. Does not include top and bottom contacts.
                           Should be N-1 regions, default is no special interfaces at all
        :param contacts: List of Tag name of the contacts for this device. Default is one at top and one at bottom
        :param environment: Local environment of the measurement. Default is 25C, no B field, 50% humidity, STP
        :param kwargs: Other kwargs for base class
        """
        super().__init__(**kwargs)
        self.regions = regions
        self.mesh = mesh
        if interfaces is None:
            interfaces = [None for _ in regions]
            interfaces.pop() # Should be 1 less interface than the number of materials.
        self.interfaces = interfaces
        if contacts is None:
            contacts = self.regions[0].mesh_data[0]['tag'], self.regions[-1].mesh_data[-1]['tag']
        self.contacts = contacts
        self.environment = environment
        for region in regions:
            logger.debug(f"Now adding region {region.name} to the new mesh for this device")
            for mesh_info in region.mesh_data:
                self.mesh.add_line(**mesh_info)
                logger.info(f"Added mesh: {mesh_info['tag']}")
            add_1d_region(mesh=self.mesh.name, material=region.material.name, region=region.name,
                          tag1=region.mesh_data[0]['tag'], tag2=region.mesh_data[-1]['tag'])
        for contact in contacts:
            logger.debug(f"Adding 1D contacts for {contact}")
            add_1d_contact(mesh=self.mesh.name, name=contact, tag=contact, material='metal')
            # TODO support contacts not ideal metals! Maybe make Contact an object? thats silly perhaps. Dict?
        self.mesh.finalize()
        logger.info(f"Creating device {self.name} associated with mesh {self.mesh.name}")
        create_device(mesh=self.mesh.name, device=self.name)

    def sweep_iv(self, start_v=-1, stop_v=1, step_v=None, count=20):
        if step_v is None:
            v_sweep = np.linspace(start_v, stop_v, count)
        else:
            v_sweep = np.arange(start_v, stop_v, step_v)
        current_sweeep = []
        for bias_volt in v_sweep:
            self.set_top_bias(bias_volt)
            current_sweeep.append(self.get_total_current())

    def set_bias(self, tag, bias_volt):
        set_parameter(device=self.name, name=f"{tag}_bias", value=float(bias_volt))

    def set_top_bias(self, bias_volt):
        tag = self.contacts[0]
        self.set_bias(tag, bias_volt)

    def set_bottom_bias(self, bias_volt):
        tag = self.contacts[1]
        self.set_bias(tag, bias_volt)

    @property
    def total_current(self):
        return self.get_total_current()

    def get_total_current(self):
        """
        Gets the total current, electron+hole from the two contacts.
        :return:
        """
        current = 0
        for contact in self.contacts:
            e_current = get_contact_current(device=self.name, contact=contact, equation=ECE_NAME)
            h_current = get_contact_current(device=self.name, contact=contact, equation=HCE_NAME)
            current += e_current + h_current
        return current
