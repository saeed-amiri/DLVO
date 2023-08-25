"""
This script initializes and displays a 2D system of particles. 
It sets their initial positions and velocities, visualizes their
distribution, and logs system-related messages.
"""

import random
from dataclasses import dataclass
import numpy as np
import scipy.stats

import logger
from read_param import ReadParam
from colors_text import TextColor as bcolors
from display_structures import DisplaySystem


# Custom exception for particle overlap
class ParticleOverlapError(Exception):
    """Raised when particles overlap after maximum placement attempts."""


@dataclass
class Particle:
    """Represents a particle in the system."""
    position: tuple[float, float]
    velocity: tuple[float, float]


class TwoDSystem:
    """prepare a 2d system"""

    def __init__(self,
                 params: dict[str, float],
                 log: logger.logging.Logger
                 ) -> None:
        self.info_msg: str = 'Messages from TwoDSystem:\n'
        self.particles = self.generate_particles(params, log)
        self.display_system_state(params)
        self.write_log_msg(log)

    def generate_initial_velocities(self,
                                    params: dict[str, float]
                                    ) -> list[tuple[float, float]]:
        """
        Set initial velocity for the particles using the
        Maxwell-Boltzmann distribution.
        """
        scale_factor = (
            2 * params['boltzmann_constant'] * params['temperature'] /
            params['particle_mass'])**0.5

        initial_velocities: list[tuple[float, float]] = []
        all_speeds = [scipy.stats.maxwell.rvs(scale=scale_factor)
                      for _ in range(int(params['num_particles']))]

        max_speed = max(all_speeds)
        normalized_speeds = [speed / max_speed for speed in all_speeds]

        for speed in normalized_speeds:
            theta = random.uniform(0, 2 * np.pi)
            vel_x = speed * np.cos(theta) * params['initial_velocities']
            vel_y = speed * np.sin(theta) * params['initial_velocities']
            initial_velocities.append((vel_x, vel_y))

        return initial_velocities

    def generate_particles(self,
                           params,
                           log: logger.logging.Logger
                           ) -> list["Particle"]:
        """initial particles"""
        initial_velocities: list[tuple[float, float]] = \
            self.generate_initial_velocities(params)
        particles: list["Particle"] = []
        occupied_positions: list[tuple[float, float]] = []
        max_try: int = 100  # Maximumn number to try to find a place
        for i in range(int(params['num_particles'])):
            try_i = 0
            while True:
                if try_i < max_try:
                    # Ensure that the particle is fully inside the 2D system
                    # considering its radius and not overlapping
                    pos_x = random.uniform(
                        params['particle_radius'],
                        params['width'] - params['particle_radius'])
                    pos_y = random.uniform(
                        params['particle_radius'],
                        params['height'] - params['particle_radius'])
                    overlap_flag: bool = False
                    for existing_pos in occupied_positions:
                        distance = \
                            ((pos_x - existing_pos[0])**2 +
                             (pos_y - existing_pos[1])**2)**0.5
                        if distance < 2 * params['particle_radius']:
                            overlap_flag = True
                            break
                    if not overlap_flag:
                        occupied_positions.append((pos_x, pos_y))
                        break
                    try_i += 1
                else:
                    err_msg: str = ('The numbers of try excceded the max '
                                    f'limit: {max_try}.\n\tReduce the numbers'
                                    ' of particles or increase the size of '
                                    'the system!\n')
                    log.error(err_msg)
                    raise ParticleOverlapError(f'\n\t{bcolors.FAIL}{err_msg}'
                                               f'{bcolors.ENDC}')

            particle = Particle((pos_x, pos_y), initial_velocities[i])
            particles.append(particle)
        return particles

    def display_system_state(self,
                             params: dict[str, float]
                             ) -> None:
        """to disply the initial structure"""
        DisplaySystem(params, self.particles)

    def write_log_msg(self,
                      log: logger.logging.Logger  # Name of the output file
                      ) -> None:
        """writing and logging messages from methods"""
        log.info(self.info_msg)
        print(f'{bcolors.OKBLUE}{self.__module__}:\n'
              f'\t{self.info_msg}\n{bcolors.ENDC}')


if __name__ == '__main__':
    LOG = logger.setup_logger(log_name='system_initiat.log')
    parameter = ReadParam(log=LOG)
    system = TwoDSystem(parameter.parameters_dict, log=LOG)
