"""
Initializing a 2d system
"""
import random

import logger
from read_param import ReadParam
from colors_text import TextColor as bcolors


class Particle:
    """set the position and velocity"""
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity


class TwoDSystem:
    """prepare a 2d system"""

    def __init__(self,
                 params: dict[str, float],
                 log: logger.logging.Logger
                 ) -> None:
        self.info_msg: str = 'Messages from TwoDSystem:\n'
        self.particles = self.initialize_particles(params)
        self.write_log_msg(log)

    def set_velocities(self,
                       params: dict[str, float]
                       ) -> list[tuple[float, float]]:
        """set initial velocity for the particles"""
        initial_velocities: list[tuple[float, float]] = []
        for _ in range(int(params['num_particles'])):
            vel_x = random.uniform(0.1, params['initial_velocities'])
            vel_y = random.uniform(0.1, params['initial_velocities'])
            initial_velocities.append((vel_x, vel_y))
        return initial_velocities

    def initialize_particles(self,
                             params) -> list["Particle"]:
        """initial particles"""
        initial_velocities: list[tuple[float, float]] = \
            self.set_velocities(params)
        particles = []
        for i in range(int(params['num_particles'])):
            # Ensure that the particle is fully inside the 2D system
            # considering its radius
            pos_x = \
                random.uniform(params['particle_radius'],
                               params['width'] - params['particle_radius'])
            pos_y = \
                random.uniform(params['particle_radius'],
                               params['height'] - params['particle_radius'])

            particle = Particle((pos_x, pos_y), initial_velocities[i])
            particles.append(particle)

        return particles

    def display(self):
        """to disply the initial structure"""
        for particle in self.particles:
            print(f"Particle -> Position: {particle.position}, "
                  f"Velocity: {particle.velocity}")

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
    system.display()
