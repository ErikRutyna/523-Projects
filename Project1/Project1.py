import numpy as np
from scipy import sparse as sp
from scipy.sparse import linalg
import math as mp
import matplotlib.pyplot as plt


def main():
    N = [0, 1, 2, 5]
    Psi = []
    xField = []
    yField = []
    uinfField = []
    pressCoeff = []
    tangentVel = []
    dragCoeff = []

    for i in range(len(N)): 
        psiTemp = streamSolver(N[i])
        xFieldTemp, yFieldTemp, uinfFieldTemp = velocityFieldGenerator(psiTemp, N[i])
        pTemp, tTemp, CDtemp = pressureCoeff(xFieldTemp, yFieldTemp, N[i])

        # Pull them back out
        Psi.append(psiTemp)
        xField.append(xFieldTemp)
        yField.append(yFieldTemp)
        uinfField.append(uinfFieldTemp)
        pressCoeff.append(pTemp)
        tangentVel.append(tTemp)
        dragCoeff.append(CDtemp)

    # Map all the contour plots
    contourGenerator(Psi[1], xField[1], yField[1], uinfField[1], N[1])

    # Tangential Mapping
    x = [[np.linspace(0, 1, 17)], [np.linspace(0, 1, 33)], [np.linspace(0, 1, 65)], [np.linspace(0, 1, 513)]]
    figure = plt.figure(figsize=(8, 6))
    plt.plot(np.transpose(x[0]), tangentVel[0])
    plt.plot(np.transpose(x[1]), tangentVel[1])
    plt.plot(np.transpose(x[2]), tangentVel[2])
    # plt.plot(np.transpose(x[3]), tangentVel[3])
    plt.legend(['r=0', 'r=1', 'r=2'])
    plt.xlabel('x', fontsize=13)
    plt.ylabel(r'Outer Wall Tangential Velocity - $u_w$', fontsize=13)
    # plt.title(r'Tangential Wall Velocity ($u_w$)', fontsize=13)
    plt.savefig('TangentialVelocity.png')

    # Velocity Profile
    h25 = 0.2 + 0.1*(1-mp.cos(mp.pi*0.25))
    h75 = 0.2 + 0.1*(1-mp.cos(mp.pi*0.75))

    # Build vectors for quiver plot
    xVec25 = np.zeros(17)
    uVec25 = np.linspace(0, h25, 17)
    yVec25 = xField[2][16]
    vVec25 = np.zeros(17)

    figure.clear()
    plt.plot(xField[2][16], np.linspace(0, h25, 17))
    # plt.quiver(xVec25, yVec25, uVec25, vVec25, width=0.005, headwidth=2, angles='xy', scale_units='xy', scale=1, color='r')
    plt.ylabel('Height', fontsize=13)
    plt.xlabel('x-Velocity', fontsize=13)
    plt.savefig('u_y_x=0.25.png')

    xVec75 = np.linspace(0, h75, 17)
    uVec75 = np.ones(17) * h75
    yVec75 = xField[2][48]
    vVec75 = np.zeros(17)
    figure.clear()
    plt.plot(xField[2][48], np.linspace(0, h75, 17))
    # plt.quiver(xVec75, yVec75, uVec75, vVec75, width=0.005, headwidth=2, angles='xy', scale_units='xy', scale=1, color='r')
    plt.ylabel('Height', fontsize=13)
    plt.xlabel('x-Velocity', fontsize=13)
    plt.savefig('u_y_x=0.75.png')

    # Convergence Study
    figure.clear()
    x = (([mp.sqrt(1/16*1/4), mp.sqrt(1/32*1/8), mp.sqrt(1/64*1/16)]))
    y = ((abs(dragCoeff[0:3] - dragCoeff[3])))
    plt.loglog(x, (y))
    plt.xlabel(r'h = $\sqrt{\Delta x \Delta y}$', fontsize=13)
    plt.ylabel(r'$|c_d - c_{d_{exact}}|$', fontsize=13)
    plt.autoscale()
    # plt.title(r'Error in Drag Coefficient', fontsize=13)
    plt.savefig('ConvergenceStudy.png', bbox_inches="tight")
    print(np.polyfit(np.log2(x), np.log2(y), 1))


def streamSolver(ref_index=0):
    """Sets up and solves the stream function system Ax = b

    Parameters
    ----------
    ref_index : int, optional
        Refinement index, increases the size of grid by powers of 2, by default 0

    Returns
    -------
    PsiM : np.array
        2-D array that forms the stream function matrix at each point in the reference domain
    """
    N = 2**(2+ref_index)
    Nxi = (4*N + 1)  # Number of grid points in xi-direction
    Neta = (N + 1)  # Number of grid points in eta-direction
    Ntot = Nxi * Neta  # Number of grid points total
    NNZ = 2*Nxi + 2*(Neta-2) + (Nxi - 1)*(Neta - 1)*9  # Number of non-zero grid points

    deta = 1 / N
    dxi = 1 / N / 4

    # For construction of the A-matrix
    rows = np.zeros(NNZ, dtype=int)
    cols = np.zeros(NNZ, dtype=int)
    data = np.zeros(NNZ, dtype=float)
    
    # Psi-value vector
    F = np.zeros(Ntot)

    iNZ = 0  # First index for non-zero value

    for iGlob in range(Ntot):
        # Grab the coordinate spaces
        Eta = mp.floor(iGlob/Nxi) / (Neta-1)
        Xi = (iGlob % Nxi) / (Nxi-1)

        # Global coordinates
        X = Xi
        h = 0.2 + 0.1 * (1 - mp.cos(mp.pi * X))
        Y = h * Eta
        
        # All of the various derivatives
        hprime = mp.pi / 10 * mp.sin(mp.pi * X)
        hdprime = mp.pi ** 2 / 10 * mp.cos(mp.pi * X)

        Etaprimex = -Y * hprime / (h ** 2)
        Etadprimexx = Y * (2 * (hprime ** 2) - hdprime * h) / (h ** 3)

        Etaprimey = 1 / h

        Xiprimex = 1

        # Check for the BC's first
        # Top floor, psi = 1
        if iGlob in range(Ntot-Nxi, Ntot):
            F[iGlob] = 1
            rows[iNZ] = cols[iNZ] = iGlob
            data[iNZ] = 1
            iNZ += 1
            continue
        # Bottom floor, psi = 0
        elif iGlob in range(Nxi-1):
            F[iGlob] = 0
            rows[iNZ] = cols[iNZ] = iGlob
            data[iNZ] = 1
            iNZ += 1
            continue
        # Left wall, psi = varies linearly
        elif iGlob % Nxi == 0:
            F[iGlob] = mp.floor(iGlob / Nxi) / N
            rows[iNZ] = cols[iNZ] = iGlob
            data[iNZ] = 1
            iNZ += 1
            continue
        # Right wall, psi = varies linearly
        elif iGlob % Nxi == (Nxi - 1):
            F[iGlob] = mp.floor(iGlob / Nxi) / N
            rows[iNZ] = cols[iNZ] = iGlob
            data[iNZ] = 1
            iNZ += 1
            continue
        # If it wasn't a BC, then it has to be an interior node

        # Assemble in row for the second derivative outside of that
        iL = iGlob - 1; iR = iGlob + 1 # Left & Right
        iD = iGlob - Nxi; iU = iGlob + Nxi # Up and Down
        iLD = iD - 1; iRD = iD + 1 # Bottom diagonals
        iLU = iU - 1; iRU = iU + 1 # Top diagonals

        # Slice for our data/col/row vectors
        localRange = range(iNZ+0, iNZ+9)

        # Locations for the A-matrix
        rows[localRange] = iGlob
        cols[localRange] = np.array([iL, iR, iU, iD, iLU, iRU, iLD, iRD, iGlob])

        # Coefficients on the Psi-terms for the A-matrix
        L = R = Xiprimex ** 2 / dxi ** 2
        U = (Etaprimex ** 2 / deta ** 2) + (Etaprimey ** 2 / deta ** 2) + (Etadprimexx / (2*deta))
        D = (Etaprimex ** 2 / deta ** 2) + (Etaprimey ** 2 / deta ** 2) - (Etadprimexx / (2*deta))
        LD = RU = 2 * Xiprimex * Etaprimex / (4 * deta * dxi)
        LU = RD = -2 * Xiprimex * Etaprimex / (4 * deta * dxi)
        Center = -2 * Xiprimex ** 2 / dxi ** 2 + (-2 * Etaprimex ** 2 / deta ** 2) + (-2 * Etaprimey ** 2 / deta ** 2)
        data[localRange] = np.array([L, R, U, D, LU, RU, LD, RD, Center])

        # data[localRange] = np.array([Xiprimex**2 / dxi**2, Xiprimex**2 / dxi**2, # Left & Right
        #     (Etaprimex ** 2 / deta ** 2) + (Etadprimexx / (2*deta)) + (Etaprimey / deta ** 2), # Up
        #     (Etaprimex ** 2 / deta ** 2) - (Etadprimexx / (2*deta)) + (Etaprimey / deta ** 2), # Down
        #     (-2*Etaprimex*Xiprimex) / (4 * deta * dxi), (2*Etaprimex*Xiprimex) / (4 * deta * dxi),# LU, RU
        #     (2 * Etaprimex * Xiprimex) / (4 * deta * dxi), (-2*Etaprimex*Xiprimex) / (4 * deta * dxi), # LD, RD
        #     (-2*Xiprimex ** 2 / dxi ** 2) + (-2 * Etaprimex ** 2 / deta ** 2) + (-2 * Etaprimey ** 2 / deta ** 2)]) # Center

        iNZ += 9

    # Build the sparse matrix up from what we have
    A = sp.csr_matrix((data, (rows, cols)), shape=(Ntot, Ntot))

    Psi = sp.linalg.spsolve(A, F)
    return Psi


def velocityFieldGenerator(Psi, ref_index=0):
    """Generates the x-velocity field and y-velocity field for a given
    stream function matrix

    Parameters
    ----------
    Psi : np.array float
        Stream function values at every grid point

    ref_index : int, optional
        Refinement index, increases the size of grid by powers of 2, by default 0

    Returns
    -------
    xVelField : np.array of float
        Velocity field in the x-direction

    yVelField : np.array of float
        Velocity field in the y-direction

    uinfVelField : np.array of float
        Velocity field magnitude
    """
    N = 2**(2+ref_index)
    Nxi = (4*N + 1) # Number of grid points in xi-direction
    Neta = (N + 1) # Number of grid points in eta-direction

    deta = 1 / N
    dxi = 1 / N / 4

    # Our 3 velocity fields
    xVelField = np.zeros((Nxi, Neta))
    yVelField = np.zeros((Nxi, Neta))
    uinfVelField = np.zeros((Nxi, Neta))

    Psi = np.reshape(Psi, (Nxi, Neta), order='F')

    # x-Velocity Field (only one we actually care about)
    for iY in range(Neta):
        for iX in range(Nxi):
            # Grab the coordinate spaces
            Eta = iY * deta
            Xi = iX * dxi

            # Global coordinates
            X = Xi
            h = 0.2 + 0.1 * (1 - mp.cos(mp.pi * X))
            Y = h * Eta
            
            # All of the various derivatives
            hprime = mp.pi / 10 * mp.sin(mp.pi * X)

            Etaprimex = -Y * hprime / (h ** 2)

            Etaprimey = 1 / h

            # Corner BC's first

            # Bottom floor
            if iY == 0:
                xVelField[iX][iY] = secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimey
                continue

            # Top top floor
            elif iY == Neta-1:
                xVelField[iX][iY] = secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimey
                continue
            
            # Anything else is an interior node
            xVelField[iX][iY] = centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimey

    # y-Velocity Field (I guess)
    for iY in range(Neta):
        for iX in range(Nxi):
            # Grab the coordinate spaces
            Eta = iY * deta
            Xi = iX * dxi

            # Global coordinates
            X = Xi
            h = 0.2 + 0.1 * (1 - mp.cos(mp.pi * X))
            Y = h * Eta
            
            # All of the various derivatives
            hprime = mp.pi / 10 * mp.sin(mp.pi * X)

            Etaprimex = -Y * hprime / (h ** 2)

            Etaprimey = 1 / h

            # Wall BC's first
            
            # Left Wall
            if iX == 0:
                if iY == 0:
                    yVelField[iX][iY] = -(secOrderRight(Psi[iX][iY], Psi[iX+1][iY], Psi[iX+2][iY], dxi) +
                        secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimex )
                    continue
                elif iY == Neta-1:
                    yVelField[iX][iY] = -(secOrderRight(Psi[iX][iY], Psi[iX+1][iY], Psi[iX+2][iY], dxi) +
                        secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimex )                    
                    continue 
                else:
                    yVelField[iX][iY] = -(secOrderRight(Psi[iX][iY], Psi[iX+1][iY], Psi[iX+2][iY], dxi) +
                        centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimex)
                    continue

            # Right Wall
            elif iX == Nxi-1:
                if iY == 0:
                    yVelField[iX][iY] = -(secOrderLeft(Psi[iX][iY], Psi[iX-1][iY], Psi[iX-2][iY], dxi) +
                        secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimex )
                    continue
                elif iY == Neta-1:
                    yVelField[iX][iY] = -(secOrderLeft(Psi[iX][iY], Psi[iX-1][iY], Psi[iX-2][iY], dxi) +
                        secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimex )                    
                    continue  
                else:
                    yVelField[iX][iY] = -(secOrderLeft(Psi[iX][iY], Psi[iX-1][iY], Psi[iX-2][iY], dxi) +
                        centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimex)
                    continue

            # Top Wall
            elif iY == Neta-1:
                yVelField[iX][iY] = -(centralDiff(Psi[iX+1][iY], Psi[iX-1][iY], dxi) +
                    secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimex )
                continue

            # Bottom Wall
            elif iY == 0:
                yVelField[iX][iY] = -(centralDiff(Psi[iX+1][iY], Psi[iX-1][iY], dxi) +
                    secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimex )
                continue

            # Everything else must be interior node
            yVelField[iX][iY] = -(centralDiff(Psi[iX+1][iY], Psi[iX-1][iY], dxi) +
                centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimex)


    # Then loop again to make the magnitude velocity field.
    for iY in range(Neta):
        for iX in range(Nxi):
            uinfVelField[iX][iY] = mp.sqrt(xVelField[iX][iY] ** 2 + yVelField[iX][iY] ** 2)

    return xVelField, yVelField, uinfVelField


def contourGenerator(Psi, xField, yField, uinfField, ref_index=0):
    """Plots the associated Contours for whichever set is given.

    Extended Summary
    ----------------
    Uses pyplot in order to plot the points used in the grid. The size
    of the grid is (4N + 1)x(N+1) in size, where N = 2^(2+ref_index). This 
    grid is then saved to a nearby out-put directory for access. 

    Parameters
    ----------
    Psi : np.array of float
        Values of the stream function at X/Y coordinates
    xField : np.array of float
        Values of the x-velocity component at X/Y coordinates
    xField : np.array of float
        Values of the y-velocity component at X/Y coordinates
    xField : np.array of float
        Values of Uinf at X/Y coordinates
    ref_index : int, optional
        Refinement index, increases the size of grid by powers of 2, by default 0
    """
    N = 2**(2+ref_index)
    Nxi = (4*N + 1) # Number of grid points in xi-direction
    Neta = (N + 1) # Number of grid points in eta-direction

    # Reference Grid
    Xi, Eta = np.meshgrid(np.linspace(0, 1, Nxi), np.linspace(0, 1, Neta))

    # Global Grid
    X = Xi
    Y = (0.2 + 0.1 * (1 - np.cos(np.pi * Xi))) * Eta


    # Plot all the coordinate pairs
    figure = plt.figure(figsize=(8, 6))
    plt.scatter(X, Y)
    plt.gca().set_aspect('equal')
    plt.xlabel('x', fontsize=13)
    plt.ylabel('y', fontsize=13)
    plt.title('Global Space Visualization', fontsize=13)
    plt.savefig('Global.png')

    figure.clear()
    plt.scatter(Xi, Eta)
    plt.gca().set_aspect('equal')
    plt.xlabel(r'$\xi$', fontsize=13)
    plt.ylabel(r'$\eta$', fontsize=13)
    plt.title('Reference Space Visualization', fontsize=13)
    plt.savefig('Ref.png')

    # Contour Mapping
    Psi = np.reshape(Psi, (Nxi, Neta), order='F')
    xField = np.reshape(xField, (Nxi, Neta), order='F')
    yField = np.reshape(yField, (Nxi, Neta), order='F')
    infField = np.reshape(uinfField, (Nxi, Neta), order='F')

    figure.clear()
    plt.contourf(X, Y, np.transpose(Psi), 9)
    plt.xlabel('x', fontsize=13)
    plt.ylabel('y', fontsize=13)
    # plt.title('$\Psi$ Contour Mapping', fontsize=13)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('StreamFunctionContour.png')

    figure.clear()
    plt.contourf(X, Y, np.transpose(xField), 20)
    plt.xlabel('x', fontsize=13)
    plt.ylabel('y', fontsize=13)
    # plt.title('x-Velocity Contour Mapping', fontsize=13)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('x-VelocityContour.png')

    figure.clear()
    plt.contourf(X, Y, np.transpose(yField), 20)
    plt.xlabel('x', fontsize=13)
    plt.ylabel('y', fontsize=13)
    # plt.title('y-Velocity Contour Mapping', fontsize=13)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('y-VelocityContour.png')

    figure.clear()
    plt.contourf(X, Y, np.transpose(infField), 20)
    plt.xlabel('x', fontsize=13)
    plt.ylabel('y', fontsize=13)
    # plt.title('$U_{\infty}$ Contour Mapping', fontsize=13)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('UinfContour.png')


def pressureCoeff(xField, yField, ref_index=0):
    """Generates the pressure coefficient for tangential along the top wall of the diffuser.

    Parameters
    ----------
    xField : np.array type float
        The x-velocities at each point in the domain

    ref_index : int, optional
        Refinement index, increases the size of grid by powers of 2, by default 0
    Returns
    -------
    c_p : np.array type float
        Array of pressure coefficients along the top wall
    
    uTangent : np.array type float
        Tangential velocity along the top wall

    CD : float
        Drag coefficient on the diffuser wall
    """
    N = 2**(2+ref_index)
    Nxi = (4*N + 1) # Number of grid points in xi-direction
    Neta = (N + 1) # Number of grid points in eta-direction

    deta = 1 / N
    dxi = 1 / N / 4

    c_p = np.zeros((Nxi))
    uTangent = np.zeros((Nxi))
    norm = np.zeros((Nxi))
    X = np.zeros((Nxi))
    CD = 0

    for iX in range(Nxi):
        # Reference space -> Global space
        Xi = iX * dxi
        X[iX] = Xi

        # Slope at each point forming a [1, h'] tangent vector
        hprime = mp.pi / 10 * mp.sin(mp.pi * X[iX])

        # Normalize the tanget vector into being a unit tangent
        magnitude = mp.sqrt(1 + hprime ** 2)
        tangentVec = np.array([1 / magnitude, hprime/magnitude])

        # Velocity vector [u, v]
        velocityVec = np.array([xField[iX][Neta-1], yField[iX][Neta-1]])

        # Tangent via dot producting
        uTangent[iX] = np.dot(tangentVec, velocityVec)

        # Formula for pressure coefficient
        c_p[iX] = 1 - (uTangent[iX] / xField[0][0]) ** 2

        # Values at various x-points for computing C_d
        h = 0.2 + 0.1*(1-np.cos(np.pi*X[iX]))
        y = h

        # Nabla * Eta
        delEtax = -y * hprime / h ** 2

        # Normalizing the x-component of the Nable * Eta vector
        norm[iX] = delEtax / mp.sqrt(delEtax ** 2 + (1/h) ** 2)

    # Integrand for trapz/by hand integration
    integrand = c_p * norm

    # This doesn't match anything, c_d ~ 0.09 for any N
    # for iX in range(Nxi-1):
    #     CD += (integrand[iX] + integrand[iX+1]) / 2 * dxi
    CD = np.trapz(integrand, x=X)

    return c_p, uTangent, CD


def wallDragCoeff(pressure, ref_index=0):
    """Integrates over the pressure coefficient to find the drag along diffuser wall.

    Parameters
    ----------
    pressure : np.array type float
        Pressure coefficients (c_p) along the top wall of the diffuser
    ref_index : int, optional
        Refinement index, by default 0

    Returns
    -------
    C_D : float
        Total drag coefficient for the diffuser
    """
    N = 2**(2+ref_index)
    Nxi = (4*N + 1) # Number of grid points in xi-direction
    Neta = (N + 1) # Number of grid points in eta-direction

    deta = 1 / N
    dxi = 1 / N / 4

    CD = 0

    # Trapezoidal Integration
    for iX in range(Nxi-1):
        # Grab the coordinate spaces
        Xi = iX * dxi

        # Global coordinates
        X = Xi
        h = 0.2 + 0.1 * (1 - mp.cos(mp.pi * X))
        Y = h
        
        # All of the various derivatives
        hprime = mp.pi / 10 * mp.sin(mp.pi * X)

        Etaprimex = -Y * hprime / (h ** 2)
        Etaprimey = 1 / h
        magnitude = (Etaprimex ** 2 + Etaprimey) ** 0.5
        Nx = Etaprimex / magnitude

        CD += (pressure[iX] + pressure[iX+1]) / 2 * dxi * Nx

    return CD


def centralDiff(uplus, uminus, delta):
    """Returns the central difference at a point

    Parameters
    ----------
    uplus : float
        Next index's value
    uminus : float
        Previous index's value
    delta : float
        Spacing between points

    Returns
    -------
        du/d(delta) at that point
    """
    return (uplus - uminus) / (2*delta)


def secOrderLeft(ought, uminus1, uminus2, delta):
    """Returns the first derivative using a left-sided schema

    Parameters
    ----------
    ought : float
        Point we're taking the derivative at
    uminus1 : float
        Previous index's value
    uminus2 : float
        2 Previous index's value
    delta : float
        Spacing between points

    Returns
    -------
        du/d(delta) at that point
    """
    return (3*ought - 4*uminus1 + uminus2) / (2*delta)


def secOrderRight(ought, uplus1, uplus2, delta):
    """Returns the first derivative using a right-sided schema

    Parameters
    ----------
    ought : float
        Point we're taking the derivative at
    uplus1 : float
        Next index's value
    uplus2 : float
        2 Next index's value
    delta : float
        Spacing between points

    Returns
    -------
        du/d(delta) at that point
    """
    return (-3*ought + 4*uplus1 - uplus2) / (2*delta)


if __name__ == "__main__":
    main()