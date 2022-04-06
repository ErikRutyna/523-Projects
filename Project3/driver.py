import timemarching as tm
import time
# This is the main driver file that runs everything for the project.


def main():
    t = time.time()
    N = 4  # [32, 32, 64, 128, 128, 512, 1024]
    Re = 10  # [100, 400, 1000, 3200, 5000, 7500, 10000]
    uWall = 1
    tol = 10**-5
    beta = 0.42

    tm.TimeMarch(N, Re, uWall, beta, tol)

    print("The simulation converged in {0} seconds.".format(time.time() - t))

    return

if __name__ == "__main__":
    main()
