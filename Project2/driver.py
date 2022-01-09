
import helper as hlp
import flux as flux
import plotmesh as pltmsh
import readgri as rdmsh

import numpy as np
import matplotlib.pyplot as plt



def main():
    CFL = 0.9
    tMax = 0.5

    meshFile = 'tank0.gri'
    # meshFile = 'tank1.gri'
    # meshFile = 'test.gri'

    # Read in the mesh
    meshInfo = rdmsh.readgri(meshFile)

    # X & Y center coordinates for each cell
    [Xc, Yc] = hlp.centroid(meshInfo['E'], meshInfo['V'])

    # Initial condition for each cell based on X-Y coordinate of the centroid
    [h0, V0] = hlp.initialCondition(Xc, Yc)

    # Calculate the initial state vector
    state0 = np.zeros((len(h0), 3))
    state0[:, 0] = h0

    # Uncomment to run free-stream state, h = 1, u = v = 0
    state0 = hlp.initialFreestream(Xc, Yc)

    # Plot the initial condition
    hlp.plotConditionState(meshInfo, h0, 'ICPlotFine.png')

    # Run the CFD for the initial condition
    stateN, stateHalf, wall, pipe1, pipe2, pipe3, dt, residuals = \
        hlp.timeIntegration(meshInfo, state0, CFL, tMax*5, tMax/2)

    # Save the state vectors and pipe forces over time
    np.savetxt('finalStateFine.csv', stateN, delimiter=',', fmt='%f')
    np.savetxt('halfwayStateFine.csv', stateHalf, delimiter=',', fmt='%f')
    np.savetxt('wallforceFine.csv', wall, delimiter=',', fmt='%f')
    np.savetxt('pipe1forceFine.csv', pipe1, delimiter=',', fmt='%f')
    np.savetxt('pipe2forceFine.csv', pipe2, delimiter=',', fmt='%f')
    np.savetxt('pipe3forceFine.csv', pipe3, delimiter=',', fmt='%f')
    np.savetxt('dtFine.csv', dt, delimiter=',', fmt='%f')
    np.savetxt('residualsFine.csv', residuals, delimiter=',', fmt='%f')

    hlp.plotConditionState(meshInfo, stateN[:, 0], 'FCPlotFine.png')
    # print('Holder Print')

    # Last minute plot for task 3
    f = plt.figure()
    plt.plot(range(len(dt)), residuals, color='blue')
    plt.xlabel('Timestep Number')
    plt.ylabel('L_2 Residual Norm')
    plt.savefig('freestream.png')
    plt.close(f)
    return


if __name__ == "__main__":
    main()