"""
For the DLVO model, the interaction energy between two colloidal
particles is given by:

U(r) = UvdW(r) + Uel(r)U(r)

Where:

    UvdW is the van der Waals interaction energy.
    Uel is the electrostatic interaction energy.

The corresponding forces are the negative gradient of these energy
profiles with respect to distance.
"""

import typing
import numpy as np

if typing.TYPE_CHECKING:
    from system_initialization import Particle


class DLVO:
    """
    A class representing the DLVO interactions.

    The DLVO theory describes the interaction between charged particles
    in a colloidal system.
    The total interaction potential is the sum of the van der Waals
    attraction and the electrostatic repulsion between the charged
    particles.

    Attributes:
        hamaker(float): Hamaker constant representing van der Waals
                        interactions.
        epsilon (float): Dielectric constant of the medium.
        psi (float): Surface potential, representing charge on particle's
                     surface.
        kappa (float): Inverse Debye length, related to ionic strength
                       of medium.
        a (float): Particle radius.

    Methods:
        compute_forces(particles: list) -> np.ndarray:
            Computes the pairwise DLVO forces on each particle.
    """
    def __init__(self,
                 params: dict[str, float],  # Parameters
                 particles: list["Particle"]  # Particles
                 ) -> None:
        self.hamaker = params['hamaker_constant']
        self.epsilon = params['dielectric_constant']
        self.psi = params['zeta_potential']
        self.kappa = params['ionic_strength']**0.5  # Inverse Debye length
        self.radius = params['particle_radius']
        self.forces = self.compute_forces(particles)

    def _vdw_force(self,
                   sep_distance: float  # Separation between particles
                   ) -> float:
        """Van der Waals force for a given separation sep_distance."""
        term: float = (2 - 2 * self.radius / sep_distance +
                       3 * self.radius**2 / sep_distance**2)
        return -self.hamaker / (6 * sep_distance**2) * term

    def _electrostatic_force(self,
                             sep_distance: float  # Separation between particle
                             ) -> float:
        """Electrostatic force for a given separation sep_distance."""
        exp_term: np.float64 = \
            np.exp(-self.kappa * (sep_distance - 2 * self.radius))
        return -2 * np.pi * self.epsilon * self.radius * self.psi**2 *\
            (self.kappa * exp_term) / (1 + exp_term)

    def total_force(self,
                    sep_distance: float  # Separation between particles
                    ) -> float:
        """Total DLVO force for a given separation sep_distance."""
        return self._vdw_force(sep_distance) + \
            self._electrostatic_force(sep_distance)

    def compute_forces(self,
                       particles: list["Particle"]  # Particles
                       ):
        """
    Compute the pairwise DLVO forces for a list of particles.

    Algorithm:
        1. Extract particle positions and form a 2D array of shape
           (N, N, 2), where N is the number of particles.
        2. Compute the pairwise position differences, delta, between
           all particles.
        3. Calculate the Euclidean distances between all particle pairs
           using the delta array.
        4. Subtract twice the particle radius from the distances to get
           the distances between particle surfaces.
        5. Create a mask to avoid computing self-interaction (force of
           a particle on itself).
        6. Using the surface distances and the mask, calculate the DLVO
           force magnitudes between particle pairs.
        7. Determine the direction of the force for each particle pair
           using the delta and distance arrays.
        8. Multiply the force magnitudes by their respective directions
           to obtain force vectors.
        9. Sum the force vectors for each particle to get the net force
           due to all other particles.

    Args:
        particles (list): List of particles with positions, where each
                           particle is expected to have a 'position'
                           attribute.

    Returns:
        np.ndarray: A 2D numpy array of shape (N, 2), where N is the
                    number of particles.
                    Each row represents the net DLVO force (in x and y
                    directions) on a particle due to all other
                    particles.
    """
        positions = np.array([particle.position for particle in particles])
        delta: np.ndarray = \
            positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
        distances: np.ndarray = np.linalg.norm(delta, axis=-1)
        surface_distances: np.ndarray = distances - 2 * self.radius

        # Mask to avoid self-interaction
        mask: np.ndarray = np.identity(len(particles)) == 0

        # Calculate force magnitudes
        force_magnitudes: np.ndarray = \
            np.vectorize(self.total_force)(surface_distances) * mask

        # Determine force directions & calculate force vectors
        with np.errstate(divide='ignore', invalid='ignore'):
            force_directions = delta / distances[:, :, np.newaxis]
            force_directions[~np.isfinite(force_directions)] = 0
        forces: np.ndarray = np.sum(force_magnitudes[:, :, np.newaxis] *
                                    force_directions, axis=1)
        return forces

def test_dlvo_class() -> None:
    """ Define test class"""

    @dataclass
    class TestParticle:
        """define particles"""
        position: tuple[float, float]
        velocity: tuple[float, float]

    # Initialize DLVO
    params: dict[str, float] = {
        'hamaker_constant': 1e-20,
        'dielectric_constant': 80.8,
        'zeta_potential': -25e-3,  # Surface potential
        'ionic_strength': 1,
        'particle_radius': 1
    }

    # Two particles located at different positions
    particle1 = TestParticle(position=(1, 1), velocity=(0, 0))
    particle2 = TestParticle(position=(2, 1), velocity=(0, 0))
    particle3 = TestParticle(position=(2, 2), velocity=(0, 0))

    # Compute DLVO forces between the two particles
    dlvo = DLVO(params, [particle1, particle2, particle3])
    print(dlvo.forces)


if __name__ == "__main__":
    from dataclasses import dataclass
    test_dlvo_class()
