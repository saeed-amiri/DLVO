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
    from sysyem_initialization import Particle


class DLVO:
    """calculate forces on DLVO model"""
    def __init__(self,
                 params: dict[str, float],  # Parameters
                 particles: list["Particle"]  # Particles
                 ) -> None:
        self.hamaker = params['hamaker_constant']
        self.epsilon = params['dielectric_constant']
        self.psi = params['zeta_potential']
        self.kappa = params['ionic_strength']**0.5  # Inverse Debye length
        self.radius = params['particle_radius']
        self.compute_forces(particles)

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
        """Compute the pairwise DLVO forces on each particle."""
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
