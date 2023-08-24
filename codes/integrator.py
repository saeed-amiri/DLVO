import numpy as np
import force_copmutation
from system_initialization import Particle, DisplaySystem


class ParticleIntegrator:
    """Class to handle the integration of the equations of motion."""

    def __init__(self,
                 params: dict[str, float],
                 particles: "list[Particle]",
                 forces: np.ndarray
                 ) -> None:
        self.delta_t = params['time_step']
        self.mass = params['particle_mass']
        self.update_system(params, particles, forces)

    def update_system(self,
                      params: dict[str, float],
                      particles: list["Particle"],
                      forces: np.ndarray
                      ) -> None:
        """Update system of particles using Velocity Verlet integration."""
        for t_step in range(int(params['total_steps'])):
            hlf_vels, new_pos = \
                self._compute_half_velocities_and_positions(
                    particles, forces, params)
            new_vels = self._compute_new_velocities(params, new_pos, hlf_vels)
            new_structure = self._update_particle_structure(new_pos, new_vels)
            DisplaySystem(params, new_structure, out_label=f'frame{t_step}')
            dlvo = force_copmutation.DLVO(params, new_structure)
            particles = new_structure
            forces = dlvo.forces.copy()

    def _compute_half_velocities_and_positions(self,
                                               particles: list["Particle"],
                                               forces: np.ndarray,
                                               params: dict[str, float]
                                               ) -> tuple[list[np.ndarray],
                                                          list["Particle"]
                                                          ]:
        new_pos: list["Particle"] = []
        hlf_vels: list[np.ndarray] = []
        for item, particle in enumerate(particles):
            v_half_dt, new_position = \
                self.velocity_verlet(particle, forces[item], params)
            hlf_vels.append(v_half_dt)
            new_pos.append(
                Particle((new_position[0], new_position[1]), (0, 0)))
        return hlf_vels, new_pos

    def _compute_new_velocities(self,
                                params: dict[str, float],
                                new_pos: list["Particle"],
                                hlf_vels: list[np.ndarray]
                                ) -> list[np.ndarray]:
        dlvo = force_copmutation.DLVO(params, particles=new_pos)
        return self._update_velocity(dlvo.forces, hlf_vels)

    def _update_particle_structure(self,
                                   new_pos: list["Particle"],
                                   new_vels: list[np.ndarray]
                                   ) -> list["Particle"]:
        new_structure = []
        for pos, vel in zip(new_pos, new_vels):
            new_structure.append(Particle(pos.position, (vel[0], vel[1])))
        return new_structure

    def _update_velocity(self,
                         forces: np.ndarray,
                         hlf_vels: list[np.ndarray]
                         ) -> list[np.ndarray]:
        """update velocity based on the updated forces"""
        new_vels: list[np.ndarray] = []
        for item, vel in enumerate(hlf_vels):
            new_vels.append(
                vel +
                0.5 * self.delta_t *
                np.array(forces[item]) / self.mass)
        return new_vels

    def velocity_verlet(self,
                        particle: "Particle",
                        force: np.ndarray,
                        params: dict[str, float]
                        ) -> tuple[np.ndarray, list[float]]:
        """Velocity Verlet integrator."""
        v_half_dt = (
            np.array(particle.velocity) +
            0.5 * self.delta_t * np.array(force) / self.mass
        )

        new_position = self._compute_new_position(particle, v_half_dt, params)
        return v_half_dt, new_position.tolist()

    def _compute_new_position(self,
                              particle: "Particle",
                              v_half_dt: np.ndarray,
                              params: dict[str, float]
                              ) -> np.ndarray:
        """
        Compute the new position of a particle considering boundary conditions.
        """
        new_position = np.array(particle.position) + self.delta_t * v_half_dt
        for i, axis in enumerate(['width', 'height']):
            while (new_position[i] - params['particle_radius'] < 0 or
                   new_position[i] + params['particle_radius'] > params[axis]):
                if new_position[i] - params['particle_radius'] < 0:
                    new_position[i] += params[axis]
                if new_position[i] + params['particle_radius'] > params[axis]:
                    new_position[i] -= params[axis]
        return new_position
