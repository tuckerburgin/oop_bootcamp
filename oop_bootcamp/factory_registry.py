# Example of a factory method for building interface/adapter classes

import abc

_toolkits = {}

# Abstract class/interface
class TrajectoryAdapter(abc.ABC):
    @abc.abstractmethod
    def compute_center_of_mass(self):
        pass

    @abc.abstractmethod
    def compute_radius_of_gyration(self):
        pass

def register(toolkit_name, toolkit_class):
    if not issubclass(toolkit_class, TrajectoryAdapter):
        raise ValueError('{0} is not a TrajectoryAdapter'.format(toolkit_class))
    _toolkits[toolkit_name] = toolkit_class     # add to _toolkits dictionary

def trajectory_factory(trajectory_toolkit, **kwargs):
    if trajectory_toolkit not in _toolkits.keys():
        raise ValueError('Toolkit not found: ' + trajectory_toolkit)

    return _toolkits[trajectory_toolkit](kwargs['filename'])
