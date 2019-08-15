# Example of programing "to a concrete" or "to a specific". This is the "bad" way of doing things...

import sys
import factory

toolkit = sys.argv[1]

myinterface = factory.trajectory_factory(toolkit, filename='protein.pdb')

print(myinterface.compute_center_of_mass())
print(myinterface.compute_radius_of_gyration())
