# Example of interface class design pattern for the same problem as in script.py

import mdtraj as md
import MDAnalysis as mda
import numpy as np
import factory_registry

# Adapter class for MDTraj
class AdaptMDTraj(factory_registry.TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = md.load_pdb(filename)

    def compute_center_of_mass(self):
        return 10 * md.compute_center_of_mass(self.trajectory)  # factor of 10 converts nm to Å

    def compute_radius_of_gyration(self):
        return 10 * md.compute_rg(self.trajectory)

# Adapter class for MDAnalysis
class AdaptMDAnalysis(factory_registry.TrajectoryAdapter):
    def __init__(self, filename):
        self.universe = mda.Universe('protein.pdb')

    def compute_center_of_mass(self):
        mass_by_frame = np.ndarray(shape=(len(self.universe.trajectory), 3))
        for ts in self.universe.trajectory:
            mass_by_frame[ts.frame] = self.universe.atoms.center_of_mass(compound='segments')
        return mass_by_frame

    def compute_radius_of_gyration(self):
        rg_by_frame = np.ndarray(shape=(len(self.universe.trajectory)))
        for ts in self.universe.trajectory:
            rg_by_frame[ts.frame] = self.universe.atoms.radius_of_gyration(compound='segments')
        return rg_by_frame


factory_registry.register('mdtraj', AdaptMDTraj)
factory_registry.register('mdanalysis', AdaptMDAnalysis)
