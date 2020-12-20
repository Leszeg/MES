import csv
from multiprocessing import Process
from multiprocessing.managers import BaseManager

from numpy import zeros

from MES.Data import global_data
from MES.FEM_Grid import FEM_Grid

"""
__init.py__ - call FEM_Grid object for startup parameters and plot grid (parallel)
During the simulation, the matrix recalculation takes place in parallel

"""


def made(g):
    """
    Function run simulation

    """

    temps = zeros((global_data.N_B * global_data.N_H, 1))
    d = {}
    for i in range(0, int(global_data.st) + int(global_data.sst), int(global_data.sst)):
        if i == 0:
            temps = g.get_temps(True, temps)[0]
            g.calculate_matrix()
            g.to_file()
        else:
            g.calculate_matrix()
        g.matrix_aggregation()
        d[f'Time {i}'] = [min(temps), max(temps)]
        temps = g.solve_ode(temps)
    with open("Temperatures.csv", "w") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(d.keys())
        writer.writerows(zip(*d.values()))


def run():
    """
    Function calls calculate grid and plot printing parallel
    """
    BaseManager.register('mc', FEM_Grid)
    manager = BaseManager()
    manager.start()
    g = manager.mc(True, 0)
    p1 = Process(target=g.plot_grid)
    p2 = Process(target=made, args=[g])
    p1.start()
    p2.start()
    p1.join()
    p2.join()
