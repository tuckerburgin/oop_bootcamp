# Example of programing "to a concrete" or "to a specific". This is the "bad" way of doing things...

import sys
from interfaces import registry_interface   # this shouldn't need to be specific, dunno what the right approach is...
# Could do: for file in interfaces_directory: import file
import factory_registry

toolkit = sys.argv[1]

myinterface = factory_registry.trajectory_factory(toolkit, filename='protein.pdb')

print(myinterface.compute_center_of_mass())
print(myinterface.compute_radius_of_gyration())
