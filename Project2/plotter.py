import numpy as np
from numpy.lib.npyio import load
import helper as hlp
import matplotlib.pyplot as plt
import plotmesh as pltmsh
import readgri as rdmsh
import pathlib
import os

def loadFile(directory, filename):
    """Loads the data requested in the filename from the directory

    Parameters
    ----------
    directory : string
        Directory the file is in
    filename : string
        The file to be loaded

    Returns
    -------
    loadedData : Usually np.array
        The data loaded from the CSV
    """
    fileToLoad = directory + '\\' + filename
    loadedData = np.genfromtxt(fileToLoad, delimiter=',', dtype=None)

    return loadedData

def plotForces(Wall, P1, P2, P3, T, DIR, fname):
    """Plots the forces for the wall and pipes in the given direction.

    Parameters
    ----------
    Mesh : dictionary
        The mesh
    Wall : np.array
        Forces along wall
    P1 : np.array
        Forces along Pipe 1
    P2 : np.array
        Forces along Pipe 2
    P3 : np.array
        Forces along Pipe 3
    T : np.array
        Time
    DIR : np.array
        X/Y direction
    fname : np.array
        Filename to be saved as
    """
    # Figure out if plotting X/Y
    index = 2
    if DIR == 'X':
        index = 0
    else: index = 1

    f = plt.figure()
    f.tight_layout()
    # plt.plot(T, Wall[:, index], color='purple')
    plt.plot(T, P1[:, index], color='red')
    plt.plot(T, P2[:, index], color='blue')
    plt.plot(T, P3[:, index], color='green')
    plt.legend(['Pipe 1','Pipe 2','Pipe 3'], loc='upper left')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')
    plt.savefig(fname)
    plt.close(f)


    return

def main():
    # Makes all the plots from the data that is already saved. Saves plots
    # in the respective data set's directory
    wd = pathlib.Path().resolve()
    stringPath = str(wd)

    # Load in all of our coarse data files
    stringCoarse = stringPath + r'\Coarse Results'

    coarseStateFinal = loadFile(stringCoarse, 'finalStateCoarse.csv')
    coarseStateHalf = loadFile(stringCoarse, 'halfwayStateCoarse.csv')
    coarseForceWall = loadFile(stringCoarse, 'wallforceCoarse.csv')
    coarseForcePipe1 = loadFile(stringCoarse, 'pipe1forceCoarse.csv')
    coarseForcePipe2 = loadFile(stringCoarse, 'pipe2forceCoarse.csv')
    coarseForcePipe3 = loadFile(stringCoarse, 'pipe3forceCoarse.csv')
    coarseResiduals = loadFile(stringCoarse, 'residualsCoarse.csv')
    coarseTimesteps = loadFile(stringCoarse, 'dtCoarse.csv')

    # Load in all of the fine data files
    stringFine = stringPath + r'\Fine Results'

    fineStateFinal = loadFile(stringFine, 'finalStateFine.csv')
    fineStateHalf = loadFile(stringFine, 'halfwayStateFine.csv')
    fineForceWall = loadFile(stringFine, 'wallforceFine.csv')
    fineForcePipe1 = loadFile(stringFine, 'pipe1forceFine.csv')
    fineForcePipe2 = loadFile(stringFine, 'pipe2forceFine.csv')
    fineForcePipe3 = loadFile(stringFine, 'pipe3forceFine.csv')
    fineResiduals = loadFile(stringFine, 'residualsFine.csv')
    fineTimesteps = loadFile(stringFine, 'dtFine.csv')

    # Mesh Infos
    meshFine = rdmsh.readgri('tank1.gri')
    meshCoarse = rdmsh.readgri('tank0.gri')

    # Change directory to this folder for saving plots
    os.chdir(stringCoarse)

    # State vector components at t=tMax
    hlp.plotConditionState(meshCoarse, coarseStateFinal[:, 0], 'coarseHeightFinal.png')
    hlp.plotConditionState(meshCoarse, coarseStateFinal[:, 1] / coarseStateFinal[:, 0], 'coarseXVelocityFinal.png')
    hlp.plotConditionState(meshCoarse, coarseStateFinal[:, 2] / coarseStateFinal[:, 0], 'coarseYVelocityFinal.png')

    # State vector components at t=tMax/2
    hlp.plotConditionState(meshCoarse, coarseStateHalf[:, 0], 'coarseHeightHalf.png')
    hlp.plotConditionState(meshCoarse, coarseStateHalf[:, 1] / coarseStateHalf[:, 0], 'coarseXVelocityHalf.png')
    hlp.plotConditionState(meshCoarse, coarseStateHalf[:, 2] / coarseStateHalf[:, 0], 'coarseYVelocityHalf.png')

    plotForces(coarseForceWall, coarseForcePipe1, coarseForcePipe2,\
        coarseForcePipe3, coarseTimesteps, 'X', 'coarseForcesX.png')
    plotForces(coarseForceWall, coarseForcePipe1, coarseForcePipe2,\
        coarseForcePipe3, coarseTimesteps, 'Y', 'coarseForcesY.png')


    # Change directory to this folder for saving plots
    os.chdir(stringFine)

    # State vector components at t=tMax
    hlp.plotConditionState(meshFine, fineStateFinal[:, 0], 'fineHeightFinal.png')
    hlp.plotConditionState(meshFine, fineStateFinal[:, 1] / fineStateFinal[:, 0], 'fineXVelocityFinal.png')
    hlp.plotConditionState(meshFine, fineStateFinal[:, 2] / fineStateFinal[:, 0], 'fineYVelocityFinal.png')

    # State vector components at t=tMax/2
    hlp.plotConditionState(meshFine, fineStateHalf[:, 0], 'fineHeightHalf.png')
    hlp.plotConditionState(meshFine, fineStateHalf[:, 1] / fineStateHalf[:, 0], 'fineXVelocityHalf.png')
    hlp.plotConditionState(meshFine, fineStateHalf[:, 2] / fineStateHalf[:, 0], 'fineYVelocityHalf.png')

    plotForces(fineForceWall, fineForcePipe1, fineForcePipe2,\
        fineForcePipe3, fineTimesteps, 'X', 'fineForcesX.png')
    plotForces(fineForceWall, fineForcePipe1, fineForcePipe2,\
        fineForcePipe3, fineTimesteps, 'Y', 'fineForcesY.png')

    return

if __name__ == "__main__":
    main()
