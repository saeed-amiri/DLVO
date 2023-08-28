"""class of method for rescaling velocity"""

import typing
import numpy as np

if typing.TYPE_CHECKING:
    from system_initialization import Particle


class Rescaling:
    """
    Different methods for rescaling the velocity of the partciles
    """
    def __init__(self,
                 particles: list["Particle"],
                 desired_temperature: float,
                 mass: float
                 ) -> None:
        self.mass = mass
        self.particles = particles
        self.desired_temperature = desired_temperature

    def scale_velocities_to_temperature(self) -> list["Particle"]:
        """
        A simple velocity rescaling method. It's one of the most
        straightforward ways to adjust the temperature of a molecular
        system. Essentially, the velocities of all the particles are
        uniformly scaled based on the ratio of the desired kinetic
        energy (related to the desired temperature) to the current
        kinetic energy.
        """
        updated_partciles: list["Particle"] = self.particles.copy()
        # Constants
        k_boltzmann: float = 1.380649e-23  # Boltzmann constant in J/K
        # Compute the current kinetic energy
        kinetik_e: float = sum(0.5 * self.mass * (v[0]**2 + v[1]**2)
                               for v in [particle.velocity
                               for particle in self.particles]
                               )
        # Calculate the scaling factor
        particle_nr: int = len(self.particles)
        scaling_fac: float = \
            np.sqrt((
                3 * (particle_nr - 1) * k_boltzmann * self.desired_temperature
                ) / (2 * kinetik_e))
        # Scale the velocities
        for particle in updated_partciles:
            particle.velocity = \
                (particle.velocity[0] * scaling_fac,
                 particle.velocity[1] * scaling_fac)
        return updated_partciles
