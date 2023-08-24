import numpy as np
import force_copmutation
from system_initialization import Particle, DisplaySystem


class ParticleIntegrator:
    """
    Class to handle the integration of the equations of motion.

    Attributes:
        dt (float): Time step for integration.
        mass (float): Mass of the particles.
    """
    def __init__(self, params, particles, forces):
        self.dt = params['time_step']
        self.mass = params['particle_mass']
        self.update_system(params, particles, forces)

    def update_system(self, params, particles, forces):
        """Update system of particles using Velocity Verlet integration."""

        for t_step in range(int(params['total_steps'])):
            new_pos = []
            new_vels = []
            hlf_vels = []
            new_structre = []
            for particle, force in zip(particles, forces):
                v_half_dt, new_position = \
                    self.velocity_verlet(particle, force, params)
                hlf_vels.append(v_half_dt)
                # Update particle position
                new_pos.append(Particle(new_position, (0, 0)))

            dlvo = force_copmutation.DLVO(params, particles=new_pos)
            for vel, force in zip(hlf_vels, dlvo.forces):
                # Update velocity with the second half of Verlet integration
                new_vels.append(
                    vel + 0.5 * self.dt * np.array(force) / self.mass)

            for pos, vel in zip(new_pos, new_vels):
                new_structre.append(Particle(pos.position, vel))
            DisplaySystem(params, new_structre, out_label=f'frame{t_step}')
            dlvo = force_copmutation.DLVO(params, new_structre)
            particles = new_structre
            forces = dlvo.forces

    def velocity_verlet(self, particle, force, params):
        """Velocity Verlet integrator."""
        # Update velocity by half a time step
        v_half_dt = (
            np.array(particle.velocity) +
            0.5 * self.dt * np.array(force) / self.mass
        )

        # Update particle position
        new_position = np.array(particle.position) + self.dt * v_half_dt
        for i, axis in enumerate(['width', 'height']):
            if new_position[i] - params['particle_radius'] < 0:
                while new_position[i] - params['particle_radius'] < 0:
                    new_position[i] += params[axis]
            if new_position[i] + params['particle_radius'] > params[axis]:
                while (
                    new_position[i] + params['particle_radius'] > params[axis]
                      ):
                    new_position[i] -= params[axis]

        # Return half-updated velocity and new position
        return v_half_dt, new_position.tolist()
