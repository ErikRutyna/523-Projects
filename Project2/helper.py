import matplotlib.pyplot as plt
import numpy as np
import math as math
import flux as flux
import readgri as rdmsh
import flux
import copy as copy

def centroid(cells, vertices):
    """Returns an array of X and Y coordinates for each cell.

    :param cells: NP array of vertex coordinate locations
    :param vertices: NP array of elements with their vertices
    :return: X, Y - NP arrays that contain X-Y coordinate pairs of "centroid" of the cell
    """
    X = np.zeros((len(cells)))
    Y = np.zeros((len(cells)))

    for K in range(len(cells)):
        # X & Y coordinates for each cell's vertices
        x = vertices[cells[K]][:, 0]
        y = vertices[cells[K]][:, 1]

        # Forumla for centroid of a triangle:
        # https://byjus.com/maths/centroid-of-a-triangle/
        X[K] = np.sum(x) / len(x)
        Y[K] = np.sum(y) / len(y)

    return X, Y

def initialCondition(X, Y):
    """Initializes the heights for the fuel in the tank based on the
    2-D height function. Velocities are set to 0

    :param X: NP Array of X-coordinates for each cell
    :param Y: NP Array of Y-coordinates for each cell
    :return: height, velocity - initial heights and velocities of the fluid for each cell
    """
    height = np.zeros((len(X)))
    velocity = np.zeros((len(X), 2))

    for K in range(len(height)):
        height[K] = 1 + 0.3 * math.exp(-50 * (X[K] - 1.3)**2 - 50 * (Y[K] - 0.9)**2)


    return height, velocity

def initialFreestream(X, Y):
    """Initializes everything to be "freestream" meaning no height and
    no velocity for anything in the fuel tank anywhere

    :param X: NP Array of X-coordinates for each cell
    :param Y: NP Array of Y-coordinates for each cell
    :return: icFreestream - initial heights and velocities at each cell
    """
    # Set everything to zero
    icFreestream = np.zeros((len(X), 3))

    # Set the heights to 1
    icFreestream[:, 0] = 1

    return icFreestream

def plotConditionState(Mesh, height, filename):
    """ Plots the given state value condition on the given mesh

    :param Mesh: Mesh from .gri file
    :param height: Height at every centroid
    :param filename: Filename to be saved as
    :return: Nothing - saves the file
    """
    # Straight up "borrowed" from plotmesh.py
    V = Mesh['V']; E = Mesh['E']; BE = Mesh['BE']
    f = plt.figure(figsize=(17.5,12))
    f.tight_layout()
    plt.tripcolor(V[:, 0], V[:, 1], triangles=E, facecolors=height, shading='flat')
    plt.set_cmap('jet')
    plt.colorbar(orientation='horizontal', pad=0.1, fraction=0.045)
    plt.tick_params(axis='both', labelsize=12)
    plt.xlabel('x-Position (m)', fontsize=28)
    plt.ylabel('y-Position (m)', fontsize=28)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(filename)
    plt.close(f)


    return

def edgePropertiesCalculator(nodeA, nodeB):
    """ Calculates the length and CCW norm out of the edge

    :param nodeA: X-Y Coordinates of node A
    :param nodeB: X-Y Coordinates of node B
    :return length: Length of the edge from A->B
    :return norm: Normal vector out of the edge in CCW fashion: [nx, ny]
    """
    # try:
    #     X = nodeA[1]
    #     Y = nodeB[1]
    # except IndexError:
    #     print(nodeA, nodeB)
    # Page 140 in KFid's notes
    length = np.linalg.norm(nodeA - nodeB)
    norm = [(nodeB[1] - nodeA[1]) / length, (nodeA[0] - nodeB[0]) / length]

    return length, norm

def areaCalculator(Mesh, cellIndex):
    """Calculates the area of the two triangular cells for the given indices.

    :param Mesh: The mesh of the problem holding all cells.
    :param cellIndex: Index of the cell
    :return A: Area of the cell
    """
    nodes = Mesh['E'][cellIndex]
    # https://en.wikipedia.org/wiki/Heron%27s_formula
    # print(cellIndex)
    a, _ = edgePropertiesCalculator(Mesh['V'][nodes[0]], Mesh['V'][nodes[1]])
    b, _ = edgePropertiesCalculator(Mesh['V'][nodes[1]], Mesh['V'][nodes[2]])
    c, _ = edgePropertiesCalculator(Mesh['V'][nodes[2]], Mesh['V'][nodes[0]])

    s = (a + b + c) / 2

    A = math.sqrt( s * (s - a) * (s - b) * (s - c) )
    return A

def sCalculator(stateVector, normal):
    """Calculates the maximum value of |s| along each edge.

    :param stateVector: Local state vector for the cell
    :param normal: Normal along the edge of the cell
    :return: s1: Maximum value of s = U +/- c where c = sqrt(g*h)
    """
    g = 9.8
    # Components of our state vector
    h = stateVector[0]
    u = stateVector[1] / stateVector[0]
    v = stateVector[2] / stateVector[0]

    # try:
    #     X = math.sqrt(g * h)
    # except ValueError:
    #     print(h)

    # Local speed
    c = math.sqrt(g * h)

    # Propagation speed
    s = max([abs(c + np.linalg.norm([u, v])), abs(c - np.linalg.norm([u, v]))])

    return s

def sCalculatorKFID(state, normalVector):
    """Re-runs KFid's |s| calculation across an edge for a given cell.

    :param state: Local state at cell i
    :param normalVector: Normal vector across edge e for cell i
    :return smax: Maximum characteristic speed along the cell
    """
    g = 9.8
    h = state[0]
    u = state[1] / h
    v = state[2] / h
    normalVelocity = np.dot([u, v], normalVector)

    c = math.sqrt(g*h)

    # eigenvalues
    l = np.zeros(3)
    l[0] = normalVelocity; l[1] = normalVelocity-c; l[2] = normalVelocity+c

    # entropy fix
    epsilon = c*.05
    for i in range(3):
        if ((l[i]<epsilon) and (l[i]>-epsilon)):
            l[i] = 0.5*(epsilon + l[i]*l[i]/epsilon)

    # absolute values of eigenvalues
    l = abs(l)

    # Get the largest eigenvalue
    smax = max(l)

    return smax

def timestepCalculator(Mesh, state, cellIndex, CFL):
    """Calculates the timestep for the given cell at its current state

    :param Mesh: Current mesh
    :param state: Current state
    :param cellIndex: Cell we're trying to find the timestep at
    :param CFL: CFL number
    :return dt: Timestep for the given cellIndex
    """
    # Local cell's area
    area = areaCalculator(Mesh, cellIndex)

    # Find the length of the 3 local edges
    edgeLengths = np.zeros((3))
    edgeNormals = np.zeros((3, 2))

    nodeA = Mesh['V'][Mesh['E'][cellIndex][0]]
    nodeB = Mesh['V'][Mesh['E'][cellIndex][1]]
    nodeC = Mesh['V'][Mesh['E'][cellIndex][2]]

    edgeLengths[0], edgeNormals[0, :] = edgePropertiesCalculator(nodeA, nodeB)
    edgeLengths[1], edgeNormals[1, :] = edgePropertiesCalculator(nodeB, nodeC)
    edgeLengths[2], edgeNormals[2, :] = edgePropertiesCalculator(nodeC, nodeA)

    # Find |s| for each side length
    s = np.zeros((3))
    s[0] = sCalculatorKFID(state[cellIndex], edgeNormals[0])
    s[1] = sCalculatorKFID(state[cellIndex], edgeNormals[1])
    s[2] = sCalculatorKFID(state[cellIndex], edgeNormals[2])

    dt = 2 * area / (np.dot(edgeLengths, s)) * CFL

    return dt

def boundaryForce(Mesh, state):
    """Returns the total force along the boundary edges at time t

    :param Mesh: Current working mesh
    :param state: State vector at the current time
    :return forces: Total force along each boundary
    """
    wallForce = np.array([0.0, 0.0])
    pipe1Force = np.array([0.0, 0.0])
    pipe2Force = np.array([0.0, 0.0])
    pipe3Force = np.array([0.0, 0.0])
    g = 9.8

    # Loop over the boundary edges and add to the forces accordingly
    for K in range(len(Mesh['BE'])):
        # Get information about the boundary edge
        nodeA = Mesh['V'][Mesh['BE'][K][0]]
        nodeB = Mesh['V'][Mesh['BE'][K][1]]
        index = Mesh['BE'][K][2]
        edgeLength, edgeNorm = edgePropertiesCalculator(nodeA, nodeB)

        # Check to see which edge it is, and add to forces accordingly
        if Mesh['BE'][K][3] == 0:
            # Wall BE
            wallForce[0] += g * state[index][0]**2 / 2 * edgeNorm[0] * edgeLength
            wallForce[1] += g * state[index][0]**2 / 2 * edgeNorm[1] * edgeLength
        elif Mesh['BE'][K][3] == 1:
            # Pipe 1 BE
            pipe1Force[0] += g * state[index][0]**2 / 2 * edgeNorm[0] * edgeLength
            pipe1Force[1] += g * state[index][0]**2 / 2 * edgeNorm[1] * edgeLength
        elif Mesh['BE'][K][3] == 2:
            # Pipe 2 BE
            pipe2Force[0] += g * state[index][0]**2 / 2 * edgeNorm[0] * edgeLength
            pipe2Force[1] += g * state[index][0]**2 / 2 * edgeNorm[1] * edgeLength
        elif Mesh['BE'][K][3] == 3:
            # Pipe 3 BE
            pipe3Force[0] += g * state[index][0]**2 / 2 * edgeNorm[0] * edgeLength
            pipe3Force[1] += g * state[index][0]**2 / 2 * edgeNorm[1] * edgeLength

    return wallForce, pipe1Force, pipe2Force, pipe3Force

def residualCalculator(Mesh, currentState, CFL, state):
    """Uses FE to march forward in time the state vector using the FVM.

    :param Mesh: Mesh of our model
    :param currentState: Initial state vector, u0
    :param CFL: The CFL number for the problem
    :return residuals: Residual vector calculator for the current state
    """
    residuals = np.zeros((len(currentState), 3))
    dt = np.zeros((len(Mesh['IE'])+len(Mesh['BE'])))

    g = 9.8

    # Loop over all the internal edges in the mesh
    for K in range(len(Mesh['IE'])):
        # [Node 1, Node 2, Left Cell Index, Right Cell Index]
        tempEdge = Mesh['IE'][K]

        # Two nodes for the edge
        nodeA = Mesh['V'][tempEdge[0]]
        nodeB = Mesh['V'][tempEdge[1]]

        # Calculate edge properties: length & norm
        edgeLength, edgeNorm = edgePropertiesCalculator(nodeA, nodeB)
        # print(tempEdge[2])

        # Areas of the two cells on the internal edge
        # A1 = areaCalculator(Mesh, tempEdge[2])
        # A2 = areaCalculator(Mesh, tempEdge[3])

        # Flux and |s| values from the two cells
        localFlux, s = flux.FluxFunction(currentState[tempEdge[2]],
                                         currentState[tempEdge[3]], edgeNorm)

        dt[K] = timestepCalculator(Mesh, state, tempEdge[2], CFL)

        # Slot the fluxes into the residuals vector
        residuals[tempEdge[2]] += localFlux * edgeLength
        residuals[tempEdge[3]] -= localFlux * edgeLength
        # print(residuals)
    # print(np.linalg.norm(residuals))
    # print('\n')
    for K in range(len(Mesh['BE'])):
        # [Node 1, Node 2, Cell Index, Boundary Name Index]
        tempEdge = Mesh['BE'][K]

        # Two nodes for the edge
        nodeA = Mesh['V'][tempEdge[0]]
        nodeB = Mesh['V'][tempEdge[1]]

        # Calculate edge properties: length & norm
        edgeLength, edgeNorm = edgePropertiesCalculator(nodeA, nodeB)
        # print(tempEdge[2])

        # Areas of the two cells on the internal edge
        # A3 = areaCalculator(Mesh, tempEdge[2])

        # Calculate the pressure terms from the BE's
        localFlux = np.zeros((3))
        # try:
        #     X = edgeNorm[0] * g * currentState[tempEdge[2]][1] ** 2 / 2
        #     Y = edgeNorm[1] * g * currentState[tempEdge[2]][1] ** 2 / 2
        # except ValueError:
        #     print(nodeA, nodeB)
        localFlux[1] = edgeNorm[0] * g * currentState[tempEdge[2]][0] ** 2 / 2
        localFlux[2] = edgeNorm[1] * g * currentState[tempEdge[2]][0] ** 2 / 2

        # s = sCalculator(state[tempEdge[2]], edgeNorm)
        
        dt[len(Mesh['IE']) + K] = timestepCalculator(Mesh, state, tempEdge[2], CFL)

        # print(residuals[tempEdge[2]])
        residuals[tempEdge[2]] += localFlux * edgeLength
        # print(residuals)
        # print('\n')
        # print(residuals)

    # print(residuals)
    # print(np.linalg.norm(residuals))
    return residuals, min(dt)

def timeIntegration(Mesh, state, CFL, maxTime, outputTime):
    """Integrates the mesh from the current state to the maximum time using forward Euler.

    :param Mesh: The mesh in place
    :param state: The state vector for all cells
    :param CFL: The CFL number for the problem
    :param maxTime: The total time to pass doing forward integration
    :param outputTime: Additional time to be output and saved
    :return finalState: The final state after doing all the time-based integration
    :return additionalStateOutput: State at the additional output time
    :return forces: Forces along boundary conditions
    """
    additionalStateOutput = state * 0
    # denominator = np.zeros((len(Mesh['E'])))
    areas = np.zeros((len(Mesh['E'])))
    # edges = np.zeros((len(Mesh['E']), 3))
    # edgeNorms = [[None for _ in range(3)] for _ in range(len(Mesh['E']))]
    # s = np.zeros((len(Mesh['E']), 3))
    wallForce = []
    pipe1Force = []
    pipe2Force = []
    pipe3Force = []
    timesteps = []
    residualNorm = []

    t = 0

    for K in range(len(Mesh['E'])):
        # Calculate the geometric properties of each cell: areas & edge lengths
        areas[K] = areaCalculator(Mesh, K)
        # print(K)
        # edges[K][0], edgeNorms[K][0] = edgePropertiesCalculator(Mesh['V'][Mesh['E']][K][0],
        #                                                         Mesh['V'][Mesh['E']][K][1])
        # edges[K][1], edgeNorms[K][1] = edgePropertiesCalculator(Mesh['V'][Mesh['E']][K][1],
        #                                                         Mesh['V'][Mesh['E']][K][2])
        # edges[K][2], edgeNorms[K][2] = edgePropertiesCalculator(Mesh['V'][Mesh['E']][K][2],
        #                                                         Mesh['V'][Mesh['E']][K][0])

    dt = 0
    # Begin integrating in time
    while t < maxTime:

        # Extra break condition in case we somehow over-march
        if t >= maxTime: break
        # Make sure we end at exactly the max time
        if (t + dt) > maxTime:
            dt = abs(maxTime - t)

        # Grab the state at a different time
        if t > 0.95 * outputTime and t < 1.05 * outputTime:
            additionalStateOutput = state

        # # Run over each cell and compute the dt
        # for K in range(len(Mesh['E'])):
        #     # Check the |s| value on each side
        #     s[K][0], s[K][1], s[K][2] = sCalculator(state[K], edgeNorms[K])
        #
        #     # Computes the denominator for the timestep formula
        #     denominator[K] = np.dot(s[K], edges[K])
        #
        #     dt[K] = 2 * areas[K] / denominator[K] * CFL

        # Calculate the force on the pipe at the current time
        wallt, pipe1t, pipe2t, pipe3t = boundaryForce(Mesh, state)

        wallForce.append(wallt)
        pipe1Force.append(pipe1t)
        pipe2Force.append(pipe2t)
        pipe3Force.append(pipe3t)

        # Our residuals and timestep
        residuals, dt = residualCalculator(Mesh, state, CFL, state)

        # Residuals and timestep outputs for plotting purposes
        timesteps.append(t)
        residualNorm.append(np.linalg.norm(residuals))

        # Some matrix-vector bullshit to make the multiplication work
        Constant = np.zeros((len(residuals), 3))
        Constant[:, 0] = dt / areas
        Constant[:, 1] = dt / areas
        Constant[:, 2] = dt / areas

        # Only uncomment this if running freestream condition
        # if np.linalg.norm(residuals) > 10**-10: print("Residual not machine precision")

        # Apply out FE Timestep
        state = state - Constant * residuals

        t += dt

    # Bring back out the final time
    finalState = state
    return finalState, additionalStateOutput, wallForce, pipe1Force, pipe2Force, pipe3Force, timesteps, residualNorm
