"""Custom Verlet scheme."""
import numpy as np
from ase.md.verlet import VelocityVerlet

class CustomVerlet(VelocityVerlet):
    """Custom Verlet scheme for ASE and FLARE.
    """
    def __init__(self, atoms, timestep=None, trajectory=None, logfile=None,
                 loginterval=1, dt=None, append_trajectory=False):

        VelocityVerlet.__init__(self, atoms, timestep, trajectory, logfile,
                                   loginterval, dt,
                                   append_trajectory=append_trajectory)
        self.old_forces = None

    def step(self, forces=None):

        atoms = self.atoms

        if forces is None:
            forces = atoms.get_forces(md=True)

        p = atoms.get_momenta()

        if self.old_forces is None:
            p = p
            #p += 0.5 * self.dt * forces
        else:
            p += 0.5 * self.dt * (forces + self.old_forces)

        atoms.set_momenta(p)

        masses = atoms.get_masses()[:, np.newaxis]
        r = atoms.get_positions()

        atoms.set_positions(r + self.dt * p / masses + 0.5 * self.dt**2 * forces / masses)

        self.old_forces = forces

        return forces
