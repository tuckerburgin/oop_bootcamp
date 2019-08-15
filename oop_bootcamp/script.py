# Example of programing "to a concrete" or "to a specific". This is the "bad" way of doing things...

import mdtraj as md
import MDAnalysis as mda
import numpy as np
import sys

toolkit = sys.argv[1]

if toolkit == 'mdtraj':         # mdtraj style
    trajectory = md.load_pdb('protein.pdb')
    print(10 * md.compute_center_of_mass(trajectory))   # factor of 10 converts nm to Å
    print(10 * md.compute_rg(trajectory))
elif toolkit == 'mdanalysis':   # MDAnalysis style
    universe = mda.Universe('protein.pdb')
    mass_by_frame = np.ndarray(shape=(len(universe.trajectory), 3))
    for ts in universe.trajectory:
        mass_by_frame[ts.frame] = universe.atoms.center_of_mass(compound='segments')
    print(mass_by_frame)
    rg_by_frame = np.ndarray(shape=(len(universe.trajectory)))
    for ts in universe.trajectory:
        rg_by_frame[ts.frame] = universe.atoms.radius_of_gyration(compound='segments')
    print(rg_by_frame)
else:
    raise AttributeError
