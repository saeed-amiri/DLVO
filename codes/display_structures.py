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
            scale_factor = 5e-3
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
