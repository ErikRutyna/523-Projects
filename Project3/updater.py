import numpy as np
import helper as hlp
from functools import partial
# The 3 functions that update WRT time:
# momentum halfstep, PPE, divergence correction


def FluxCalc(uField, vField, pField, Re):
    """Half-step updates the velocity via the momentum equation.

    :param uField: u-velocity field
    :param vField: v-velocity field
    :param pField: pressure field
    :param Re: Reynolds Number
    :param pool: multiprocessing pool
    :return: uFieldHalf, vFieldHalf - u & v fields updated to half-step
    """
    N = len(pField) - 2
    # Re = u*L*rho / nu -> nu = u*L*rho / Re -> nu = 1 / Re
    nu = 1.0 / Re

    # Get an array of all the indices that are going to be updated
    cellIndices = hlp.FluxIndicesGen(N)
    nodeIndices = hlp.HFluxIndicesGen(N)

    F = np.zeros((N+2, N+2))
    G = np.zeros((N+2, N+2))

    HX = np.zeros((N+3, N+3))
    HY = np.zeros((N+3, N+3))

    for k in range(len(cellIndices)):
        i = cellIndices[k, 0]
        j = cellIndices[k, 1]

        F[i, j] = hlp.FCalc([i, j], nu, uField)
        G[i, j] = hlp.GCalc([i, j], nu, vField)
    print("\n")
    for k in range(len(nodeIndices)):
        i = nodeIndices[k, 0]
        j = nodeIndices[k, 1]

        HX[i, j] = hlp.HXCalc([i, j], nu, uField, vField)
        HY[i, j] = hlp.HYCalc([i, j], nu, uField, vField)

    return F, G, HX, HY


def VelHalfUpdate(U, V, F, G, HX, HY, dt):
    """ Uses the fluxes to update the velocity to the half-step time node.

    :param U: x-velocity field
    :param V: y-velocity field
    :param F: x-momentum flux in x-direction
    :param G: y-momentum flux in y-direction
    :param HX: x-momentum flux in y-direction
    :param HY: y-momentum flux in x-direction
    :param dt: timestep
    :return: uField, vField velocity fields at next half-step in time
    """
    N = len(U) - 2

    # Need to get indices of the fields
    UIndex = hlp.UIndicesGen(N)
    VIndex = hlp.VIndicesGen(N)

    for k in range(len(UIndex)):
        i = UIndex[k, 0]
        j = UIndex[k, 1]

        U[i, j] = hlp.UHalf([i, j], U, F, HX, dt)

    for k in range(len(VIndex)):
        i = VIndex[k, 0]
        j = VIndex[k, 1]

        V[i, j] = hlp.VHalf([i, j], V, G, HY, dt)

    return U, V


def PPE(uField, vField, pField, dt):
    """Does RB-GS to update the pressure field until it matches the stepped
    velocity

    :param uField: x-velocity field
    :param vField: y-velocity field
    :param pField: pressure field at n
    :param dt: timestep
    :param pool: multiprocessing pool
    :return: pField - pressure field at n+1
    """
    N = len(pField) - 2
    residual = 1
    h = 1 / N

    redList = hlp.RedIndicesGen(N)
    blackList = hlp.BlackIndicesGen(N)

    PUpdatePartial = partial(hlp.PUpdate, P=pField, U=uField, V=vField, dt=dt, h=h)

    for i in range(50):
    # while (residual > 10 ** -6):
        # Grab the updated red nodes
        for k in range(len(redList)):
            i = redList[k, 0]
            j = redList[k, 1]

            pField[i, j] = PUpdatePartial([i, j])

        # Update the ghost cells for my indexing
        # Top wall
        pField[N + 1, :] = pField[N, :]
        # Bottom wall
        pField[0, :] = pField[1, :]
        # Left wall
        pField[:, 0] = pField[:, 1]
        # Right wall
        pField[:, N + 1] = pField[:, N]

        # Grab the updated black nodes
        for k in range(len(blackList)):
            i = blackList[k, 0]
            j = blackList[k, 1]

            pField[i, j] = PUpdatePartial([i, j])

        # Update the ghost cells for my indexing
        # Top wall
        pField[N + 1, :] = pField[N, :]
        # Bottom wall
        pField[0, :] = pField[1, :]
        # Left wall
        pField[:, 0] = pField[:, 1]
        # Right wall
        pField[:, N + 1] = pField[:, N]

        # Calculate our residual
        # residual = hlp.GSResidual(pField, uField, vField, dt)

    # print(pField)
    # Reduce the mean so it doesn't get really high
    #pField -= np.mean(pField[1:N+1, 1:N+1])

    return pField


def PPE2(uField, vField, pField, dt):
    """Does GS to update the pressure field until it matches the stepped
    velocity

    :param uField: x-velocity field
    :param vField: y-velocity field
    :param pField: pressure field at n
    :param dt: timestep
    :return: pField - pressure field at n+1
    """
    N = len(pField) - 2
    residual = 1
    h = 1 / N

    PUpdatePartial = partial(hlp.PUpdate, P=pField, U=uField, V=vField, dt=dt, h=h)

    for z in range(100):
        # Grab the updated red nodes
        for i in range(1, N+1):
            for j in range(1, N+1):
                pField[i, j] = PUpdatePartial([i, j])

        # Update the ghost cells for my indexing
        # Top wall
        pField[N + 1, :] = pField[N, :]
        # Bottom wall
        pField[0, :] = pField[1, :]
        # Left wall
        pField[:, 0] = pField[:, 1]
        # Right wall
        pField[:, N + 1] = pField[:, N]

    # Reduce the mean so it doesn't get really high
    #pField -= np.mean(pField[1:N+1, 1:N+1])

    return pField


def Correction(U, V, P, dt):
    """ Corrects the velocity field from n+1/2 to n+1

    :param U: x-velocity field
    :param V: y-velocity field
    :param P: pressure field
    :param dt: timestep
    :return: U, V same fields, but updated to n+1
    """

    # Fields no longer have some weird indexing shenanigans
    N = len(P) - 2

    # Find the indices of velocities that are changing
    uVelIndices = hlp.UIndicesGen(N)
    vVelIndices = hlp.VIndicesGen(N)

    # Functions for parallelization
    uLocalPartial = partial(hlp.LocalUUpdate, U=U, P=P, dt=dt)
    vLocalPartial = partial(hlp.LocalVUpdate, V=V, P=P, dt=dt)

    # Find the updated velocities and update the matrices
    for k in range(len(uVelIndices)):
        i = uVelIndices[k, 0]
        j = uVelIndices[k, 1]

        U[i, j] = uLocalPartial([i, j])

    for k in range(len(vVelIndices)):
        i = vVelIndices[k, 0]
        j = vVelIndices[k, 1]

        V[i, j] = vLocalPartial([i, j])

    return U, V



