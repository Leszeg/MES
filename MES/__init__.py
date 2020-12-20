import csv
from multiprocessing import Process, Pool

from numpy import zeros, dot
from numpy.linalg import solve

from MES.Data import global_data
from MES.FEM_Grid import FEM_Grid

"""
__init.py__ - call FEM_Grid object for startup parameters and plot grid (parallel)
During the simulation, the matrix recalculation takes place in parallel

"""

temps = zeros((global_data.N_B * global_data.N_H, 1))
g = FEM_Grid(True, temps)


def made():
    """
    Function run simulation

    """
    temps = zeros((global_data.N_B * global_data.N_H, 1))
    d = {}
    pool = Pool()
    for i in range(0, int(global_data.st) + int(global_data.sst), int(global_data.sst)):
        global g
        if i == 0:
            g.get_temps(True, temps)
            temps = g.temperature_of_nodes
            g.to_file()
            pool.map(recalculate_matrix, range(len(g.elements)))
        else:
            pool.map(recalculate_matrix, range(len(g.elements)))
        g.matrix_aggregation()
        d[f'Time {i}'] = temps[0]

        Hz = g.H_global + (g.C_global / global_data.sst)
        x = -dot(g.C_global / global_data.sst, temps.T).T
        Pz = -(g.P_global + x)
        temps = solve(Hz, Pz.T).T

        g.get_temps(False, temps)
    with open("Temperatures.csv", "w") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(d.keys())
        writer.writerows(zip(*d.values()))


def recalculate_matrix(i):
    """
        Function replaces loops and is performed in parallel
    Parameters
    ----------
    i : int
        Iterator for elements

    Returns
    -------

    """
    g.elements[i].H_matrix()
    g.elements[i].C_matrix()
    g.elements[i].boundary_condition()


def run():
    """
    Function calls calculate grid and plot printing parallel
    """
    p1 = Process(target=made)
    p2 = Process(target=g.plot_grid)

    p1.start()
    p2.start()
    p1.join()
    p2.join()
