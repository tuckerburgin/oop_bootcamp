# Example of programing "to a concrete" or "to a specific". This is the "bad" way of doing things...

import sys
import interface    # PyCharm complains about this for some reason, but it is valid

toolkit = sys.argv[1]

if toolkit == 'mdtraj':         # mdtraj style
    myinterface = interface.AdaptMDTraj('protein.pdb')
elif toolkit == 'mdanalysis':   # MDAnalysis style
    myinterface = interface.AdaptMDAnalysis('protein.pdb')
else:
    raise ValueError('analysis package variable must be mdanalysis or mdtraj')

print(myinterface.compute_center_of_mass())
print(myinterface.compute_radius_of_gyration())
