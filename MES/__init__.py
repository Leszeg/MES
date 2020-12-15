from MES.Data import print_grid_data, global_data
from MES.FEM_Grid import FEM_Grid
import numpy as np


def get_temps(g: FEM_Grid):
    x = g.temperature_of_nodes.T
    return x.T


def mult(X, Y):
    result = np.zeros((16, 1))
    # iterate through rows of X
    for i in range(len(X)):
        # iterate through columns of Y
        for j in range(len(Y[0])):
            # iterate through rows of Y
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result


temps = np.zeros((16, 1))
for i in range(0, 500, 50):
    if i == 0:
        g = FEM_Grid(True, temps)
        temps = get_temps(g)
    else:
        g = FEM_Grid(False, temps)
        temps = get_temps(g)
    Hz = g.H_global + (g.C_global / global_data.sst)
    x = -mult(g.C_global / global_data.sst, temps.T).T
    Pz = -(g.P_global + x)
    t = np.linalg.inv(Hz)
    temps = mult(t, Pz.T).T
    pass