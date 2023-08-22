"""
Reading parametes for all the simulation
System Parameters:

    System Size:
        width: Width of the 2D system.
        height: Height of the 2D system.

    Number of Particles:
        num_particles: The total number of particles in your 2D system.

    Particle Properties:
        particle_radius: Radius of the particles. If they vary in size,
        this can be an array.
        particle_mass: Mass of the particles, if considering dynamics.

    Particle Initial Conditions:
        initial_positions: Initial positions of particles, often
                           randomly assigned.
        initial_velocities: Initial velocities of particles, which
                            could be zero or based on some distribution.

DLVO Model Parameters:

    Van der Waals Interaction:
        Hamaker_constant: Represents the strength of the van der Waals
                          interactions.

    Electrostatic Interaction:
        zeta_potential: Surface potential, which indicates the
        electrostatic repulsion between particles.
        dielectric_constant: Dielectric constant of the medium.
        ionic_strength: Measure of the concentration of ions in the
        solution.
        Boltzmann_constant: Physical constant representing the relation
        between temperature and energy.
        temperature: Temperature of the system in Kelvin.

    Solvent Properties (which might affect electrostatic interactions):
        solvent_viscosity: The viscosity of the solvent, if considering
        hydrodynamic interactions or Brownian motion.
        solvent_density: Density of the solvent.

Simulation Parameters:

    Time Integration:
        time_step: The size of the simulation time step.
        total_steps: Total number of steps the simulation will run for.

    Boundary Conditions:
        Boundary type (e.g., periodic, reflective, or open).

    Thermal Fluctuations (if considered):
        friction_coefficient: If including Langevin dynamics to simulate
                              thermal fluctuations.
        random_force_magnitude: Magnitude of random forces if
                                considering Brownian motion.

    Interaction Modifiers:
        Parameters to tune or adjust interactions, like scaling factors,
        if you plan on "tuning" interactions during the simulation.

    Analysis & Visualization:
        Frequency of data logging, visualization, or analysis.
        Types of metrics or quantities to record (e.g., potential energy,
        kinetic energy, cluster sizes).

    Optimization Parameters:
        If implementing more advanced algorithms like spatial
        decomposition (e.g., cell lists) for faster computation, then
        cell size and related parameters might be needed.

"""

import json
import typing

import logger
import my_tools
from colors_text import TextColor as bcolors


class ReadParam:
    """
    Reads and stores all the parameters from the parameter file.
    """

    # Name of the parameter file is static
    FNAME: str = 'param'
    # List of the parameters that should be given to the script
    PARAMETERS: list[str] = \
        [
            "width", "height", "num_particles", "particle_radius",
            "particle_mass", "initial_positions", "initial_velocities",
            "hamaker_constant", "zeta_potential", "dielectric_constant",
            "ionic_strength", "boltzmann_constant", "temperature",
            "solvent_viscosity", "solvent_density", "time_step",
            "total_steps", "friction_coefficient", "random_force_magnitude"
        ]

    def __init__(self,
                 log: logger.logging.Logger
                 ) -> None:
        # To save the parameters read from the file
        self.parameters_dict: dict[str, typing.Any] = {}
        self.info_msg: str = 'Messages from ReadPram:\n'
        self.read_parameter_file()
        self.write_log_msg(log)

    def read_parameter_file(self) -> None:
        """
        read parameter file one by one and save it to a dictionary
        """
        # Open the file in read mode
        with open(self.FNAME, 'r', encoding='utf8') as f_r:
            for line in f_r:
                stripped_line = line.strip()
                if stripped_line.startswith("@"):
                    key, value = self._process_line(stripped_line)
                    self.parameters_dict[key] = value
        self.info_msg += '\tParametrs are:\n'
        self.info_msg += json.dumps(self.parameters_dict, indent=8)

    @staticmethod
    def _process_line(line: str) -> list[typing.Any]:
        """break the line and return the key and the value"""
        line = my_tools.drop_string(line, '@')
        parts: list[str] = line.split('=')
        return parts

    def write_log_msg(self,
                      log: logger.logging.Logger  # Name of the output file
                      ) -> None:
        """writing and logging messages from methods"""
        log.info(self.info_msg)
        print(f'{bcolors.OKBLUE}{self.__module__}:\n'
              f'\t{self.info_msg}\n{bcolors.ENDC}')


if __name__ == '__main__':
    ReadParam(log=logger.setup_logger(log_name='read_param.log'))
