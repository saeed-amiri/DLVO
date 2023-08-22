"""
Initializing a 2d system
"""
import random
import numpy as np
import scipy.stats
import matplotlib.pylab as plt

import logger
from read_param import ReadParam
from colors_text import TextColor as bcolors


class DisplaySystem:
    """
    Print a snapshot of the system in its initial conditions
    """

    transparent: bool = False  # In saving the fig

    def __init__(self,
                 params: dict[str, float],  # Parameters read from param file
                 particles: list["Particle"]  # System with their particles
                 ) -> None:
        self.plot_structure(params, particles)
        self.plot_velocity_distribution(particles)

    def plot_structure(self,
                       params: dict[str, float],  # Parameters read from param
                       particles: list["Particle"]  # System with particles
                       ) -> None:
        """plotting the structure"""

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

    def plot_velocity_distribution(self, particles: list["Particle"]) -> None:
        """Plotting the distribution of particle velocities."""

        # Extracting velocity magnitudes from particles
        velocities = \
            [((particle.velocity[0]**2 + particle.velocity[1]**2)**0.5)
             for particle in particles]

        print(velocities)
        # Create histogram
        plt.hist(velocities, bins=50, edgecolor='black', alpha=1)
        plt.title("Particle Velocity Distribution")
        plt.xlabel("Velocity")
        plt.ylabel("Number of Particles")
        plt.grid(True)
        plt.show()

    @staticmethod
    def save_close_fig(fig: plt.figure,  # The figure to save,
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
                    transparent=transparent
                    )
        plt.close(fig)


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
        self.display(params)
        self.write_log_msg(log)

    def set_velocities(self,
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

    def initialize_particles(self,
                             params) -> list["Particle"]:
        """initial particles"""
        initial_velocities: list[tuple[float, float]] = \
            self.set_velocities(params)
        particles: list["Particle"] = []
        occupied_positions: list[tuple[float, ...]] = []
        for i in range(int(params['num_particles'])):
            while True:
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

            particle = Particle((pos_x, pos_y), initial_velocities[i])
            particles.append(particle)
        return particles

    def display(self,
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
