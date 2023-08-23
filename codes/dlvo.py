"""
Modeling DLVO
"""


import logger
import read_param
import system_initialization
import force_copmutation

if __name__ == "__main__":
    LOG = logger.setup_logger(log_name='dlvo.log')
    parameter = read_param.ReadParam(log=LOG)
    system = \
        system_initialization.TwoDSystem(parameter.parameters_dict, log=LOG)
    forces = \
        force_copmutation.DLVO(parameter.parameters_dict, system.particles)
