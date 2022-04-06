import numpy as np
import matplotlib.pyplot as plt

# This file contains all the post processing functions to plot the solutions.

def main():
    N = 32
    Re = 100

    ufname = "uVelocity" + "_" + str(N) + "_" + str(Re) + ".csv"
    vfname = "vVelocity" + "_" + str(N) + "_" + str(Re) + ".csv"

    uField = np.loadtxt(ufname, delimiter=",")
    vField = np.loadtxt(vfname, delimiter=',')

    uField = uField[1:N+1, 1:N+2]
    vField = vField[1:N+2, 1:N+1]

    uGhia = np.array([0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1])
    xGhia = [0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1]
    psi = np.zeros((N+1, N))

    for i in range(N+1):
        psiTemp = 0
        for j in range(N):
            psiTemp -= vField[i, j]
            psi[i, j] += psiTemp

    levels = ([3e-3, 1.5e-3, 1.0e-3, 5.0e-4, 2.5e-4, 1.0e-4, 5.e-5, 1.0e-5, 1.e-6, 1.e-7, 1.e-8, -1e-10, -1e07, -1e-5,\
               -1e-4, -0.01, -0.03, -0.05, -0.07, -0.09, -0.1, -0.11, -0.115, -0.1175])
    X = np.linspace(0, 1, N)
    Y = np.linspace(0, 1, N+1)
    X, Y = np.meshgrid(X, Y)
    plt.figure()
    plt.contour(X, Y, psi)
    plt.show()






if __name__ == "__main__":
    main()