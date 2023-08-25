"""
Documentation: Particle System Visualization
Overview:

This module provides tools for visualizing a 2D particle system.
It mainly comprises a class, DisplaySystem, which is equipped with
methods to display the initial structure of particles and their
velocity distribution.
Dependencies:

    typing: Used for type hints and annotations.
    matplotlib.pylab: Used for plotting purposes.

Classes:

    DisplaySystem

        Purpose: Handles the visualization aspects of a 2D particle
                 system.

        Attributes:
            out_label (str): Prefix for the names of output
                             visualization files.

        Methods:

            __init__(self, params: dict[str, float],
                     particles: list["Particle"],
                     out_label: str = 'frame_i'):
            Initializes the DisplaySystem and immediately calls
            visualization methods.

            display_initial_structure(self, params: dict[str, float],
                                      particles: list["Particle"]):
            Visualizes the initial positions and velocities of the
            particles on a 2D plane.

            display_velocity_histogram(self, particles: list["Particle"]):
            Plots a histogram showcasing the distribution of particle
            velocities.

            save_close_fig(cls, fig: plt.figure, axs: plt.axes,
                           fname: str, legend=True,
                           loc: str = 'upper right', transparent=False):
            A class method to save a given figure to an output file
            and then close it.

How to Use:

    Import the necessary classes/modules.
    Create an instance of the DisplaySystem class with appropriate
    parameters and a list of Particle objects.
    Visualizations will be automatically generated upon instantiation.
"""

import typing
import matplotlib.pylab as plt

if typing.TYPE_CHECKING:
    from system_initialization import Particle


class DisplaySystem:
    """
    Handles the visualization aspects of the 2D particle system.
    """

    def __init__(self,
                 params: dict[str, float],  # Parameters read from param file
                 particles: list["Particle"],  # System with their particles
                 out_label: str = 'frame_i'  # Initiate of the output name
                 ) -> None:
        self.out_label: str = out_label
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
            if particle.velocity[0] == 0 and particle.velocity[1] == 0:
                head_width = 0
                head_length = 0
            else:
                head_width = params['particle_radius']/2
                head_length = params['particle_radius']/3
                
            scale_factor = 5e-3
            ax_i.arrow(particle.position[0],
                       particle.position[1],
                       particle.velocity[0] * scale_factor,
                       particle.velocity[1] * scale_factor,
                       head_width=head_width,
                       head_length=head_length,
                       fc='red',
                       ec='red')

        # Enhance plot appearance
        ax_i.set_aspect('equal', 'box')
        ax_i.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax_i.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
        ax_i.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

        # Display the plot
        self.save_close_fig(
            fig_i, ax_i, f'{self.out_label}_struc.png', legend=False)

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
        self.save_close_fig(
            fig_i, ax_i, f'{self.out_label}_velo_dist.png', legend=False)

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
                    transparent=transparent
                    )
        plt.close(fig)
