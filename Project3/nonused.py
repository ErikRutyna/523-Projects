import numpy as np


def CellFluxIndicesGen(N):
    """Generates an array of zero-indexed coordinate pairs of cells that
    need their fluxes computed at

    :param N: Size of mesh in 1 direction
    :return: L: array of index pairs zero-indexed that are cells where we compute fluxes
    """
    L = []

    for i in range(N):
        for j in range(N):
            L.append([i, j])

    return np.array(L)


def NodeFluxIndicesGen(N):
    """Generates an array of zero-indexed coordinate pairs of nodes that
    need their fluxes computed at

    :param N: Size of mesh in 1 direction
    :return: L: array of index pairs zero-indexed that are nodes where we compute fluxes
    """
    L = []

    for i in range(N + 1):
        for j in range(N + 1):
            # Skip bottom and top row
            if (i == 0) or (i == N):
                continue
            # Remove left and right wall
            elif (j == 0) or (j == N):
                continue
            # Anything left must be a node bordering a cell to be updated
            else:
                L.append([i, j])

    return np.array(L)


def velCalc(i, uField, vField, direction):
    """Finds the velocity at cell i for the given direction/field.

    :param i: Cell index
    :param uField: x-velocity field
    :param vField: y-velocity field
    :param direction: direction of the velocity term
    :return: returns the respective velocity term
    """
    # For u*u calculations
    if direction == "xx":
        # Index the four possible needed velocity nodes for QUICK/SMART
        u = [uField[i[0], i[1]-1], uField[i[0], i[1]], uField[i[0], i[1]+1], uField[i[0], i[1]+2]]
        # This is our averaged velocity
        q = 0.5 * (u[1] + u[2])
        # Find apply the SMART limiting
        phi = SMART(q, u)

        return q * phi

    # For v*v calculations
    elif direction == "yy":
        # Index the four possible needed velocity nodes for QUICK/SMART
        v = [vField[i[0]-2, i[1]], vField[i[0]-1, i[1]], vField[i[0], i[1]], vField[i[0]+1, i[1]]]
        # This is our averaged velocity
        q = 0.5 * (v[1] + v[2])
        # Find apply the SMART limiting
        phi = SMART(q, v)

        return q*phi
    # For u*v calculations
    elif direction == "xy":
        # Index the four possible needed velocity nodes for QUICK/SMART
        v = [vField[i[0]+1, i[1]-1], vField[i[0]+1, i[1]+0], vField[i[0]+1, i[1]+1], vField[i[0]+1, i[1]+2]]
        # This is our averaged velocity
        q = 0.5 * (uField[i[0], i[1]+1] + uField[i[0]+1, i[1]+1])
        # Find apply the SMART limiting
        phi = SMART(q, v)

        return q * phi

    # For v*u calculations
    elif direction == "yx":
        # Index the four possible needed velocity nodes for QUICK/SMART
        u = [uField[i[0]-1, i[1]+1], uField[i[0]+0, i[1]+1], uField[i[0]+1, i[1]], uField[i[0]+2, i[1]+1]]
        # This is our averaged velocity
        q = 0.5 * (vField[i[0]+1, i[1]] + vField[i[0]+1, i[1]+1])
        # Find apply the SMART limiting
        phi = SMART(q, u)

        return q * phi


def LocalVeldUdX(i, Field, h):
    """ Central difference on local u-velocity field for
    x-velocity in a given cell

    :param i: Cell index
    :param Field: x-velocity field
    :param h: length of cell
    :return: u - local x-velocity
    """
    u = (Field[i[0], i[1]+1] - Field[i[0], i[1]]) / h
    return u


def LocalVeldUdY(i, Field, h):
    """ Central difference on local u-velocity field for
    x-velocity in a given cell

    :param i: Cell index
    :param Field: x-velocity field
    :param h: length of cell
    :return: u - local x-velocity
    """
    u = (Field[i[0], i[1]+1] - Field[i[0], i[1]]) / h
    return u


def LocalVeldVdY(i, Field, h):
    """ Central difference on local v-velocity field for
    y-velocity in a given cell

    :param i: Cell index
    :param Field: y-velocity field
    :param h: length of cell
    :return: v - local y-velocity
    """
    v = (Field[i[0], i[1]+1] - Field[i[0], i[1]]) / h
    return v


def LocalVeldVdX(i, Field, h):
    """ Central difference on local v-velocity field for
    y-velocity in a given cell

    :param i: Cell index
    :param Field: y-velocity field
    :param h: length of cell
    :return: v - local y-velocity
    """
    v = (Field[i[0]+1, i[1]] - Field[i[0], i[1]]) / h
    return v


def LocalUVelUpdate(index, field, P, dt, h):
    """ Updates the local x-velocity at the index

    :param index: index of velocity to be updated
    :param field: x/y velocity field
    :param P: pressure field
    :param dt: timestep
    :param h: distance between cells
    :return: updatedVel: locally updated velocity at n+1
    """
    updatedUVel = field[index[0]+1, index[1]+1] - dt / h * (P[index[0], index[1]+1] - P[index[0], index[1]+0])

    return updatedUVel


def LocalVVelUpdate(index, field, P, dt, h):
    """ Updates the local y-velocity at the index

    :param index: index of velocity to be updated
    :param field: x/y velocity field
    :param P: pressure field
    :param dt: timestep
    :param h: distance between cells
    :return: updatedVel: locally updated velocity at n+1
    """

    updatedVVel = field[index[0]+2, index[1]+1] - dt / h * (P[index[0]+1, index[1]+0] - P[index[0]+0, index[1]+0])

    return updatedVVel

