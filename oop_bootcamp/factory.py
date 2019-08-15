# Example of a factory method for building interface/adapter classes

import interface

def trajectory_factory(trajectory_toolkit, **kwargs):
    traj_toolkits = {'mdtraj': interface.AdaptMDTraj, 'mdanalysis': interface.AdaptMDAnalysis}

    if trajectory_toolkit not in traj_toolkits.keys():
        raise ValueError('analysis package variable must be mdanalysis or mdtraj')

    return traj_toolkits[trajectory_toolkit](kwargs['filename'])
