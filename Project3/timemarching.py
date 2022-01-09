import helper as hlp
import boundary as bdy
import updater as updt
import numpy as np
import matplotlib.pyplot as plt
# This file is strictly the time marching scheme for the project.


def TimeMarch(N, Re, wallVel, beta, tol):
    """Takes in the given parameters, discretizes the domain, then marches forward
    in time in order to reach steady state for the given parameters.

    :param N: Number of square cells in each direction
    :param Re: The Reynold's Number for the flow
    :param rho: The density of the fluid
    :param wallVel: The velocity of the upper wall
    :param beta: stability factor
    :param tol: The tolerance value for checking if the solution is at steady state
    :return: Does not return anything, but saves the final state to text for post processing.
    """
    residualNorm = 1
    residualNorms = []
    residualNormX = []
    residualNormY = []
    i = 0
    t = 0
    h = 1 / N
    nu = 1 / Re
    dt = beta * min([h ** 2 / (4 * nu), (4 * nu)])

    # Setup the initial condition of the mesh
    uField, vField, pField = bdy.IC(N, wallVel)
    # uField, vField, pField = bdy.ICTest(N)
    ufname = "uVelocity" + "_" + str(N) + "_" + str(Re) + ".csv"
    vfname = "vVelocity" + "_" + str(N) + "_" + str(Re) + ".csv"
    nfname = "errorNorm" + "_" + str(N) + "_" + str(Re) + ".csv"


    # Until tolerance is met, keep marching forward in time.
    while (residualNorm > tol):
        # Print a nice little line to console to let us know how the simulation is doing
        if (i % 25) == 0:
            print("Iteration: {0} at time t = {1} s. Residual is at: {2}".format(i, t, residualNorm))

        # Exit the simulation early if we get some sort of instability that
        # causes the residuals to explode
        if residualNorm > 100:
            print("Instability detected! solution is no longer converging!")
            print("Exiting on iteration: {0} and time t = {1} seconds w/ residuals at: {2}".format(i, t, residualNorm))
            print(residualNorms)
            break

        # Update the states at the ghost nodes before doing any sort of integration
        uField, vField, pField = bdy.GhostUpdate(uField, vField, pField, wallVel)

        # Solve for the fluxes to update the velocity terms
        F, G, HX, HY = updt.FluxCalc(uField, vField, pField, Re)

        # Update the velocity terms using the fluxes and pressures
        uField, vField = updt.VelHalfUpdate(uField, vField, F, G, HX, HY, dt)

        # Solve the PPE for the fully updated pressure field at n+1
        pField = updt.PPE2(uField, vField, pField, dt)

        # Use the n+1 pressure field to correct the velocity field to n+1
        uField, vField = updt.Correction(uField, vField, pField, dt)

        # Calculate error norm
        Rx, Ry = hlp.ErrorNorm(pField, F, G, HX, HY, h)

        residualNorm = Rx + Ry
        # Update some trackers
        i += 1
        t += dt
        residualNorms.append(residualNorm)

        # stepDivergence = hlp.DivCalc(uField, vField, N)
        # Save needed data so that we don't have to run a bunch of simulations
    np.savetxt(ufname, uField, fmt='%1.5f', delimiter=',')
    np.savetxt(vfname, vField, fmt='%1.5f', delimiter=',')
    np.savetxt(nfname, residualNorms, fmt='%1.10f', delimiter=',')
    print("Simulation Complete and outputs have been saved! Please check outputs via the Post-Processing File.")
    return