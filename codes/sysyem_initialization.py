"""
This script initializes and displays a 2D system of particles. 
It sets their initial positions and velocities, visualizes their
distribution, and logs system-related messages.
"""

import random
from dataclasses import dataclass
import numpy as np
import scipy.stats
import matplotlib.pylab as plt

import logger
from read_param import ReadParam
from colors_text import TextColor as bcolors


# Custom exception for particle overlap
class ParticleOverlapError(Exception):
    """Raised when particles overlap after maximum placement attempts."""
    pass


class DisplaySystem:
    """
    Handles the visualization aspects of the 2D particle system.
    """

    TRANSPARENT: bool = False  # In saving the fig

    def __init__(self,
                 params: dict[str, float],  # Parameters read from param file
                 particles: list["Particle"]  # System with their particles
                 ) -> None:
        self.display_initial_structure(params, particles)
        self.display_velocity_histogram(particles)

    def display_initial_structure(self,
                                  params: dict[str, float],  # Parameters
                                  particles: list["Particle"]  # System atoms
                                  ) -> None:
        """
        Visualize the initial positions and velocities of the particles.
        """

        # Set up the plot with the given dimensions
        fig_i, ax_i = plt.subplots(figsize=(10, 10))
        ax_i.set_xlim(0, params['width'])
        ax_i.set_ylim(0, params['height'])
        ax_i.set_xlabel("X-axis")
        ax_i.set_ylabel("Y-axis")
        ax_i.set_title("Initial Positions of Particles in 2D System")

        # Plot each particle
        for particle in particles:
            circle = plt.Circle(particle.position, params['particle_radius'],
                                color='blue', alpha=0.6, ec='black')
            ax_i.add_patch(circle)
            # Draw an arrow representing velocity. We scale the arrow
            # size for better visualization.
            scale_factor = 0.5
            ax_i.arrow(particle.position[0],
                       particle.position[1],
                       particle.velocity[0] * scale_factor,
                       particle.velocity[1] * scale_factor,
                       head_width=params['particle_radius']/2,
                       head_length=params['particle_radius']/3,
                       fc='red',
                       ec='red')

        # Enhance plot appearance
        ax_i.set_aspect('equal', 'box')
        ax_i.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax_i.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
        ax_i.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

        # Display the plot
        self.save_close_fig(fig_i, ax_i, 'initial_pos.png', legend=False)

    def display_velocity_histogram(self,
                                   particles: list["Particle"]
                                   ) -> None:
        """Plotting the distribution of particle velocities."""
        fig_i, ax_i = plt.subplots(figsize=(10, 10))

        # Extracting velocity magnitudes from particles
        velocities = \
            [((particle.velocity[0]**2 + particle.velocity[1]**2)**0.5)
             for particle in particles]

        # Create histogram
        ax_i.hist(velocities, bins=50, edgecolor='black', alpha=1)
        ax_i.set_title("Particle Velocity Distribution")
        ax_i.set_xlabel("Velocity")
        ax_i.set_ylabel("Number of Particles")
        ax_i.grid(True)
        self.save_close_fig(fig_i, ax_i, 'velocities_dist.png', legend=False)

    @classmethod
    def save_close_fig(cls,
                       fig: plt.figure,  # The figure to save,
                       axs: plt.axes,  # Axes to plot
                       fname: str,  # Name of the output for the fig
                       legend=True,
                       loc: str = 'upper right',  # Location of the legend
                       transparent=False,
                       ) -> None:
        """
        Save the figure and close it.

        This method saves the given figure and closes it after saving.

        Args:
            fig (plt.figure): The figure to save.
            axs (plt.axes): The axes to plot.
            fname (str): Name of the output file for the figure.
            loc (str, optional): Location of the legend. Default is
            'upper right'.
        """
        if legend:
            legend = axs.legend(loc=loc)
        fig.savefig(fname,
                    dpi=300,
                    pad_inches=0.1,
                    edgecolor='auto',
                    bbox_inches='tight',
                    transparent=cls.TRANSPARENT
                    )
        plt.close(fig)


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
