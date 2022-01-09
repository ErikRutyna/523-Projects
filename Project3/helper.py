import numpy as np
from functools import partial
import boundary as bdy
# This file contains all the helper functions used for the time marching scheme.


def UIndicesGen(N):
    """ Generates a list of the cell indexes to be used with flux calculations,
    this is a list of all the cells that get updated for the half-step since
    these are the only cells where velocity changes

    :param N: Number of cells in mesh
    :return: indexList - array of indices
    """
    indexList = []
    for i in range(1, N+1):
        for j in range(2, N+1):
            indexList.append([i, j])

    return np.array(indexList)


def VIndicesGen(N):
    """ Generates a list of the cell indexes to be used with flux calculations,
    this is a list of all the cells that get updated for the half-step since
    these are the only cells where velocity changes

    :param N: Number of cells in mesh
    :return: indexList - array of indices
    """
    indexList = []
    for i in range(2, N+1):
        for j in range(1, N+1):
            indexList.append([i, j])

    return np.array(indexList)


def SMART(phi):
    """Implements the SMART limiting scheme for a given velocity and local
    velocity values along the same axis

    :param phi: A velocity list that has the -1, 0, +1 indices
    :return: phi_half: SMART-applied local velocity
    """
    phi_j = 0
    # If the phi's equate out, can't have /0, so you get the value itself
    if phi[2] == phi[0]:
        return phi[1]

    phi_hat = (phi[1] - phi[0]) / (phi[2] - phi[0])

    # Check to see which phi_hat to use
    if phi_hat <= 0 or phi_hat >= 1:
        phi_j = phi_hat

    elif 0 < phi_hat <= (1 / 6):
        phi_j = 3 * phi_hat

    elif (1 / 6) < phi_hat <= (5 / 6):
        phi_j = 3 / 8 * (2 * phi_hat + 1)

    elif (5 / 6) < phi_hat < 1:
        phi_j = 1

    #print("Phi: {0}".format(phi_j))
    # Then update the value for phi_j
    phi_j = phi[0] + (phi[2] - phi[0]) * phi_j

    return phi_j


def FluxIndicesGen(N):

    indexList = []

    for i in range(1, N+1):
        for j in range(1, N+1):
            indexList.append([i, j])

    return np.array(indexList)


def HFluxIndicesGen(N):

    indexList = []

    for i in range(1, N+2):
        for j in range(1, N+2):
            indexList.append([i, j])

    return np.array(indexList)


def FCalc(index, nu, uField):
    """Calculates the x-momentum flux in the x-direction.

    :param index: Cell index
    :param nu: kinematic viscosity
    :param uField: u-velocity field
    :return: x-momentum in the x-direction
    """
    N = len(uField) - 2
    h = 1 / N

    # Indices for finding velocities
    i = index[0]
    j = index[1]

    # Find q and then find phi according to q
    q = (uField[i, j+1] + uField[i, j]) / 2
    if q > 0:
        phi = SMART([uField[i, j-1], uField[i, j], uField[i, j+1]])
    elif q < 0:
        phi = SMART([uField[i, j+2], uField[i, j+1], uField[i, j]])
    elif q == 0:  # If q = 0, then phi doesn't matter
        phi = 0



    F = q * phi - nu * (uField[i, j+1] - uField[i, j]) / h

    return F


def GCalc(index, nu, vField):
    """Calculates the y-momentum flux in the y-direction

    :param index: Cell index
    :param nu: kinematic viscosity
    :param vField: v-velocity field
    :return: y-momentum in the y-direction
    """
    N = len(vField) - 3
    h = 1 / N

    # Indices for finding velocities
    i = index[0]
    j = index[1]

    # Find q and then find phi according to q
    q = (vField[i+1, j] + vField[i, j]) / 2
    if q > 0:
        phi = SMART([vField[i-1, j], vField[i, j], vField[i+1, j]])
    elif q < 0:
        phi = SMART([vField[i+2, j], vField[i+1, j], vField[i, j]])
    elif q == 0:  # If q = 0, then phi doesn't matter
        phi = 0

    #print("Index: {0}\t q: {1}\t Phi: {2}".format(index, q, phi))

    G = q * phi - nu * (vField[i+1, j] - vField[i, j]) / h

    return G


def HXCalc(index, nu, uField, vField):
    """Calculates the x-momentum flux in the y-direction.

    :param index: Cell index
    :param nu: kinematic viscosity
    :param uField: u-velocity field
    :param vField: v-velocity field
    :return: x-momentum in the y-direction
    """
    N = len(uField) - 2
    h = 1 / N

    # Indices for finding velocities
    i = index[0]
    j = index[1]

    # Find q and then find phi according to q
    q = (vField[i, j] + vField[i, j-1]) / 2
    if q > 0:
        phi = SMART([uField[i-2, j], uField[i-1, j], uField[i, j]])
    elif q < 0:
        phi = SMART([uField[i+1, j], uField[i, j], uField[i-1, j]])
    elif q == 0:  # If q = 0, then phi doesn't matter
        phi = 0

    HX = q * phi - nu * (uField[i, j] - uField[i-1, j]) / h

    return HX


def HYCalc(index, nu, uField, vField):
    """Calculates the y-momentum flux in the x-direction.

    :param index: Cell index
    :param nu: kinematic viscosity
    :param uField: u-velocity field
    :param vField: v-velocity field
    :return: y-momentum in the x-direction
    """
    N = len(uField) - 2
    h = 1 / N

    # Indices for finding velocities
    i = index[0]
    j = index[1]

    # Find q and then find phi according to q
    q = (uField[i, j] + uField[i-1, j]) / 2
    if q > 0:
        phi = SMART([vField[i, j-2], vField[i, j-1], vField[i, j]])
    elif q < 0:
        phi = SMART([vField[i, j+1], vField[i, j], vField[i, j-1]])
    elif q == 0:  # If q = 0, then phi doesn't matter
        phi = 0

    HY = q * phi - nu * (vField[i, j] - vField[i, j-1]) / h

    return HY


def UHalf(index, U, F, HX, dt):
    """Returns the new velocity at the index by forward Euler updating
    via the momentum equation using the fluxes for the u-velocity field

    :param index: index for velocity to be updated
    :param U: x-Velocity field
    :param F: U-momentum in U-direction
    :param HX: u-momentum in v-direction
    :param dt: timestep
    :return: uNext: updated x-velocity at the index
    """
    N = len(U) - 2
    h = 1 / N

    # Indices
    i = index[0]
    j = index[1]

    Fx = (F[i, j] - F[i, j-1]) / h
    Hy = (HX[i+1, j] - HX[i, j]) / h

    uNext = U[i, j] - dt * (Fx + Hy)

    return uNext


def VHalf(index, V, G, HY, dt):
    """Returns the new velocity at the index by forward Euler updating
    via the momentum equation using the fluxes for the v-velocity field

    :param index: index for velocity to be updated
    :param V: y-Velocity field
    :param G: V-momentum in V-direction
    :param HY: v-momentum in u-direction
    :param dt: timestep
    :return: vNext: updated y-velocity at the index
    """
    N = len(V) - 3
    h = 1 / N

    # Indices
    i = index[0]
    j = index[1]

    Gy = (G[i, j] - G[i-1, j]) / h
    Hx = (HY[i, j+1] - HY[i, j]) / h

    vNext = V[i, j] - dt * (Gy + Hx)

    return vNext


def PUpdate(index, P, U, V, dt, h):
    """Updates the red node at index using RB-GS update

    :param index: coordinate pair of node
    :param P: Pressure field
    :param U: x-Velocity field
    :param V: y-Velocity field
    :param dt: timestep
    :param h: distance between nodes
    :return: Pnext: updated local pressure value
    """
    w = 1.6

    # Indices
    i = index[0]
    j = index[1]

    pTerm = (P[i - 1, j] + P[i + 1, j] + P[i, j - 1] + P[i, j + 1])

    divVTerm = h / dt * (U[i, j + 1] - U[i, j] + V[i + 1, j] - V[i, j])

    Pnext = w / 4 * (pTerm - divVTerm) + (1 - w) * P[i, j]

    # divV = 1 / dt * (U[i, j + 1] - U[i, j] + V[i + 1, j] - V[i, j])
    #
    # sumP = (P[i - 1, j] + P[i + 1, j] + P[i, j - 1] + P[i, j + 1] - 4 * P[i, j]) / h ** 2
    #
    # error = abs(divV - sumP)
    # print("The error on this iteration of GS is: {0}".format(error))

    return Pnext


def ErrorNorm(P, F, G, HX, HY, h):
    """ Calculates the residual at the timestep

    :param P: pressure field
    :param F: x-momentum
    :param G: y-momentum
    :param HX: y-momentum in x-direction
    :param HY: x-momentum in y-direction
    :param h: spacing
    :return: L1Norm: residual from timestep
    """
    L1Normx = 0
    L1Normy = 0
    N = len(P) - 2
    for i in range(1, N+1):
        for j in range(2, N+1):
            L1Normx += abs(h * (F[i, j] - F[i, j-1]) + h * (P[i, j] - P[i, j-1]) + h * (HX[i+1, j] - HX[i, j]))
            # L1Norm += abs(h * (P[i, j] - P[i, j-1]))
            # L1Norm += abs(h * (HX[i+1, j] - HX[i, j]))

    for i in range(2, N+1):
        for j in range(1, N+1):
            L1Normy += abs(h * (G[i, j] - G[i-1, j]) + h * (P[i, j] - P[i-1, j]) + h * (HY[i, j+1] - HY[i, j]))
            # L1Norm += abs(h * (P[i, j] - P[i-1, j]))
            # L1Norm += abs(h * (HY[i, j+1] - HY[i, j]))

    return L1Normx, L1Normy


def LocalUUpdate(index, U, P, dt):
    """ Updates the singular x-velocity node at the index to N+1

    :param index: local index
    :param U: x-velocity field
    :param P: pressure field
    :param dt: timestep
    :return: uNext
    """
    N = len(P) - 2
    h = 1 / N

    # For indexing purposes
    i = index[0]
    j = index[1]

    uNext = U[i, j] - dt / h * (P[i, j] - P[i, j-1])

    return uNext


def LocalVUpdate(index, V, P, dt):
    """ Updates the singular y-velocity node at the index to N+1

    :param index: local index
    :param V: y-velocity field
    :param P: pressure field
    :param dt: timestep
    :return: vNext
    """
    N = len(P) - 2
    h = 1 / N

    # For indexing purposes
    i = index[0]
    j = index[1]

    vNext = V[i, j] - dt / h * (P[i, j] - P[i-1, j])

    return vNext


def RedIndicesGen(N):
    """List of the nodes for the red updates in RB-GS

    :param N: Size of pressure field
    :return: indexList - array of index pairs
    """
    indexList = []

    for i in range(1, N+1):
        for j in range(1, N+1):
            if i % 2 != 0:
                if j % 2 != 0:
                    indexList.append([i, j])
                else:
                    continue
            else:
                if j % 2 == 0:
                    indexList.append([i, j])
                else:
                    continue

    return np.array(indexList)


def BlackIndicesGen(N):
    """List of the nodes for the red updates in RB-GS

    :param N: Size of pressure field
    :return: indexList - array of index pairs
    """
    indexList = []

    for i in range(1, N+1):
        for j in range(1, N+1):
            if i % 2 != 0:
                if j % 2 != 0:
                    continue
                else:
                    indexList.append([i, j])
            else:
                if j % 2 == 0:
                    continue
                else:
                    indexList.append([i, j])

    return np.array(indexList)


def PError(index, P, U, V, h, dt):
    """ Calculates error at cell index for PPE solving

    param index: cell zero-indexed index
    param P: pressure field
    param U: x-velocity field
    param V: y-velocity field
    param h: size
    param dt: timestep
    return: error: RHS - LHS
    """
    i = index[0]
    j = index[1]

    pTerm = (P[i, j-1] + P[i, j+1] + P[i-1, j] + P[i+1, j])
    velTerm = (U[i, j+1] - U[i, j] + V[i+1, j] - V[i, j]) * h * dt

    error = (velTerm - pTerm + 4 * P[i, j]) / 4

    return abs(error)


def GSResidual(P, U, V, dt):
    """ Calculates the L1 error norm for a iteration of SOR RB-GS

    param P: pressure field
    param U: x-velocity field
    param V: y-velocity field
    param dt: timestep
    param pool: multiprocessing pool
    return: Error Norm
    """
    N = len(P) - 2
    h = 1 / N

    pIndices = FluxIndicesGen(N)

    ErrorPartial = partial(PError, P=P, U=U, V=V, h=h, dt=dt)

    errors = []

    for k in range(len(pIndices)):
        i = pIndices[k, 0]
        j = pIndices[k, 1]

        errors.append(ErrorPartial([i, j]) ** 2)

    return ((sum(np.array(errors))) ** 0.5) / len(errors)


def DivCalc(U, V, N):
    """ Checks if the flow field is incompressible

    """
    h = 1 / N
    compCheck = np.zeros((N+2, N+2))
    for i in range(1, N+1):
        for j in range(1, N+1):
            UR = U[i, j+1]
            UL = U[i, j]
            VT = V[i+1, j]
            VB = V[i, j]
            compCheck[i, j] = ((UR - UL) / h + (VT - VB) / h)

    return compCheck