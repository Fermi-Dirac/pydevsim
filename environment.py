from ds import set_parameter
from adaptusim import setup_logger
logger = setup_logger(__name__)

class Environment(object):
    """
    This class is the base class for more esoteric environments

    Environment holds attributes and equations to manage all external influences on the device.
    Ambient temperature (Maybe this should be allowed to change per material...)
    External magnetic field
    Light sources
    """
    def __init__(self,temp=293, B_field=0, light=None):
        self.temp = temp
        self.B_field = B_field
        self.light = light
        # for name, value in locals().items():
        #     set_parameter(name=name, value=value)