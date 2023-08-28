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
        total_ke: float = sum(0.5 * self.mass * (v[0]**2 + v[1]**2)
                              for v in [particle.velocity
                              for particle in self.particles]
                              )
        # Calculate the scaling factor
        particle_nr: int = len(self.particles)
        scaling_fac: float = \
            np.sqrt((
                3 * (particle_nr - 1) * k_boltzmann * self.desired_temperature
                ) / (2 * total_ke))
        # Scale the velocities
        for particle in updated_partciles:
            particle.velocity = \
                (particle.velocity[0] * scaling_fac,
                 particle.velocity[1] * scaling_fac)
        return updated_partciles

    def boltzmann_rescale(self) -> list["Particle"]:
        """
        Rescale velocities of the particles to achieve the desired
        temperature based on the boltzman method.

        Args:
            particles (list[Particle]): List of particle objects.
            T_desired (float): Desired temperature for the system.

        Returns:
            list[Particle]: List of particles with rescaled velocities.
        """
        updated_particles: list["Particle"] = self.particles.copy()
        # Calculate the current kinetic energy of the system
        total_ke = 0.0
        for particle in self.particles:
            total_ke += \
                0.5 * self.mass * np.dot(particle.velocity, particle.velocity)

        # Compute the current temperature
        # k_boltzmann = 1.0
        k_boltzmann: float = 1.380649e-23  # Boltzmann constant in J/K
        particle_nr = len(self.particles)
        temp_current = (2 / (3 * particle_nr * k_boltzmann)) * total_ke

        # Compute the scaling factor
        scaling_fac = np.sqrt(self.desired_temperature / temp_current)
        # Rescale the velocities of all particles
        for particle in updated_particles:
            particle.velocity = \
                tuple(scaling_fac * np.array(particle.velocity))

        return updated_particles
