import numpy as np
from scipy import sparse as sp
from scipy.sparse import linalg
import math as mp
import matplotlib.pyplot as plt

def main():
    N = 2
    Psi = streamSolver(N)
    xField, yField, infField = velocityFieldGenerator(Psi, N)
    contourPlotter(Psi, xField, yField, infField, N)

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
    Nxi = (4*N + 1) # Number of grid points in xi-direction
    Neta = (N + 1) # Number of grid points in eta-direction
    Ntot = Nxi * Neta # Number of grid points total
    NNZ = 2*Nxi + 2*(Neta-2) + (Nxi - 1)*(Neta - 1)*9 # Number of non-zero grid points

    deta = 1 / N
    dxi = 1 / N / 4

    # For construction of the A-matrix
    rows = np.zeros(NNZ, dtype=int)
    cols = np.zeros(NNZ, dtype=int)
    data = np.zeros(NNZ, dtype=float)
    
    # Psi-value vector
    F = np.zeros(Ntot)

    iNZ = 0 # First index for non-zero value

    for iGlob in range(Ntot):
        # Grab the coordinate spaces
        Eta = mp.floor(iGlob/Nxi) / Neta
        Xi = (iGlob % Nxi + 1) / Nxi

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
        elif  iGlob % Nxi == 0:
            F[iGlob] = mp.floor(iGlob / Nxi) / N
            rows[iNZ] = cols[iNZ] = iGlob
            data[iNZ] = 1
            iNZ += 1
            continue
        # Right wall, psi = varies linearly
        elif  iGlob % Nxi == (Nxi - 1):
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
    xVelField = np.zeros((Neta, Nxi))
    yVelField = np.zeros((Neta, Nxi))
    uinfVelField = np.zeros((Neta, Nxi))

    Psi = np.reshape(Psi, (Nxi, Neta), order='F')

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

            # Corner BC's
            # Bottom left
            if iY == 0 and iX == 0:
                xVelField[iY][iX] = secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimey
                yVelField[iY][iX] = -(secOrderRight(Psi[iX][iY], Psi[iX+1][iY], Psi[iX+2][iY], dxi) + \
                    secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimex)
                continue
            
            # Bottom right
            elif iY == 0 and iX == Nxi-1:
                xVelField[iY][iX] = secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimey
                yVelField[iY][iX] = -(secOrderLeft(Psi[iX][iY], Psi[iX-1][iY], Psi[iX-2][iY], dxi) +\
                    secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimex)
                continue

            # Top left
            elif iY == Neta-1 and iX == 0:
                xVelField[iY][iX] = secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimey
                yVelField[iY][iX] = -(secOrderRight(Psi[iX][iY], Psi[iX+1][iY], Psi[iX+2][iY], dxi) + \
                    secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimex)
                continue
            
            # Top right
            elif iY == Neta-1 and iX == Nxi-1:
                xVelField[iY][iX] = secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimey
                yVelField[iY][iX] = -(secOrderLeft(Psi[iX][iY], Psi[iX-1][iY], Psi[iX-2][iY], dxi) + \
                    secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimex)           
                continue

            # Top/bottom wall BC
            if iY == 0 or iY == Neta-1:
                # Corners already done
                if iX == 0  or iX == Nxi-1:
                    continue
                elif iY == Neta-1:
                    xVelField[iY][iX] = secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimey
                    yVelField[iY][iX] = -(centralDiff(Psi[iX+1][iY], Psi[iX-1][iY], dxi) +\
                        secOrderLeft(Psi[iX][iY], Psi[iX][iY-1], Psi[iX][iY-2], deta) * Etaprimex) 
                    continue
                elif iY == 0:
                    xVelField[iY][iX] = secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimey
                    yVelField[iY][iX] = -(centralDiff(Psi[iX+1][iY], Psi[iX-1][iY], dxi) +\
                        secOrderRight(Psi[iX][iY], Psi[iX][iY+1], Psi[iX][iY+2], deta) * Etaprimex)  
                    continue

            # Left/right wall BC
            if iX == 0 or iX == Nxi-1:
                # Corners already done
                if iY == 0  or iY == Neta-1:
                    continue
                elif iX == Nxi-1:
                    xVelField[iY][iX] = centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimey
                    yVelField[iY][iX] = -(secOrderLeft(Psi[iX][iY], Psi[iX-1][iY], Psi[iX-2][iY], dxi)+\
                        centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimex)
                    continue

                elif iX == 0:
                    xVelField[iY][iX] = centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimey
                    yVelField[iY][iX] = -(secOrderRight(Psi[iX][iY], Psi[iX+1][iY], Psi[iX+2][iY], dxi)+\
                        centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimex)
                    continue

            # Anything else is interior node
            xVelField[iY][iX] = centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimey
            yVelField[iY][iX] = -(centralDiff(Psi[iX+1][iY], Psi[iX-1][iY], dxi) +\
                centralDiff(Psi[iX][iY+1], Psi[iX][iY-1], deta) * Etaprimex)

    # Then loop again to make the magnitude velocity field.
    for iY in range(Neta):
        for iX in range(Nxi):
            uinfVelField[iY][iX] = mp.sqrt(xVelField[iY][iX] ** 2 + yVelField[iY][iX] ** 2)

    return xVelField, yVelField, uinfVelField

def contourPlotter(Psi, xField, yField, infField, ref_index=0):
    """Draws a picutres of the grid (Tron fans rejoice!) & associated contours.

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
    figure = plt.figure()
    plt.scatter(X, Y)
    plt.gca().set_aspect('equal')
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('Global Space Visualization', fontsize=16)
    plt.savefig('Global.png')

    figure.clear()
    plt.scatter(Xi, Eta)
    plt.gca().set_aspect('equal')
    plt.xlabel(r'$\xi$', fontsize=16)
    plt.ylabel(r'$\eta$', fontsize=16)
    plt.title('Reference Space Visualization', fontsize=16)
    plt.savefig('Ref.png')

    # Contour Mapping
    Psi = np.reshape(Psi, (Nxi, Neta), order='F')
    xField = np.reshape(xField, (Nxi, Neta), order='F')
    yField = np.reshape(yField, (Nxi, Neta), order='F')
    infField = np.reshape(infField, (Nxi, Neta), order='F')

    figure.clear()
    plt.contourf(X, Y, np.transpose(Psi), 20)
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('$\Psi$ Contour Mapping', fontsize=16)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('Contour.png')

    figure.clear()
    plt.contourf(X, Y, np.transpose(xField))
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('x-Velocity Contour Mapping', fontsize=16)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('x-Velocity.png')

    figure.clear()
    plt.contourf(X, Y, np.transpose(yField))
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('y-Velocity Contour Mapping', fontsize=16)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('y-Velocity.png')

    figure.clear()
    plt.contourf(X, Y, np.transpose(infField))
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('$U_{\infty}$ Contour Mapping', fontsize=16)
    plt.colorbar()
    # plt.plot(X,Y,'-',color='red')
    # plt.plot(X.transpose(),Y.transpose(),'-',color='red')
    plt.savefig('Uinf.png')

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