import numpy as np
# This file has all the updates to the boundary values, so cells
# like ghost cell updating and IC's



def IC(N, wallVel):
    """Generates the initial condition based on the wall velocity of the top wall.

    :param N: Number of cells in the mesh (N x N), total mesh w/ ghost is (N+1)x(N+1)
    :param wallVel: Velocity of the top wall in the positive x-direction
    :return: uField, vField, pField: Velocity and pressure fields over the domain
    INCLUDING THE GHOST CELLS
    """
    # Initialize the fields to their certain sizes - turns out all values are 0
    uField = np.zeros((N+2, N+3))
    vField = np.zeros((N+3, N+2))
    pField = np.zeros((N+2, N+2))

    return uField, vField, pField


def GhostUpdate(uField, vField, pField, uWall):
    """Updates the values at the ghost nodes.

    :param uField: u-velocity field matrix
    :param vField: v-velocity field matrix
    :param pField: pressure field matrix
    :param uWall: velocity along the top wall
    :return: uField, vField, pField - but with updated ghost nodes
    """
    N = len(pField) - 2

    # u-Field - 4 unique wall BCs

    # Apply the no-flow through BC enforcement before checking/updating ghost cells
    uField[:, 1] = 0
    uField[:, N+1] = 0

    # Left wall - copy inside
    uField[:, 0] = uField[:, 2]
    # Right Wall - copy inside
    uField[:, N+2] = uField[:, N]
    # Top wall - apply odd function
    uField[N+1, :] = 2 * uWall - uField[N, :]
    # Bottom wall - apply odd function
    uField[0, :] = 2 * 0 - uField[1, :]


    # v-Field - 4 unique wall BCs

    # Apply the no-flow through BC enforcement before checking/updating ghost cells
    vField[1, :] = 0
    vField[N+1, :] = 0

    # Top wall - copy inside
    vField[N+2, :] = vField[N, :]
    # Bottom wall - copy inside
    vField[0, :] = vField[2, :]
    # Left wall - apply odd function
    vField[:, 0] = 2 * 0 - vField[:, 1]
    # Right wall - apply odd function
    vField[:, N+1] = 2 * 0 - vField[:, N+0]

    # Pressure-Field - just copy over the nearest in-domain neighbor
    # Top wall
    pField[N+1, :] = pField[N, :]
    # Bottom wall
    pField[0, :] = pField[1, :]
    # Left wall
    pField[:, 0] = pField[:, 1]
    # Right wall
    pField[:, N+1] = pField[:, N]

    return uField, vField, pField


def ICTest(N):
    """ Makes IC but for debugging purposes

    param N: Number of interier domain cells
    :return: uField, vField, pField: Velocity and pressure fields over the domain
    INCLUDING THE GHOST CELLS
    """
    # Initialize the fields to their certain sizes - turns out all values are 0
    h = 1 / N
    uField = np.zeros((N+2, N+3))
    vField = np.zeros((N+3, N+2))
    pField = np.zeros((N+2, N+2))

    for i in range(1, N+2):
        uField[:, i] = (i-1) * h

    for i in range(1, N+2):
        vField[i, :] = -(i-1) * h

    for i in range(1, N+2):
        pField[:, i] = np.exp(i)

    return uField, vField, pField