import csv

from numpy import zeros

from MES.Data import global_data
from MES.FEM_Grid import FEM_Grid

"""
__init.py__ - call FEM_Grid object for startup parameters
"""


def made(g):
    """
    Function run simulation
    """

    temps = zeros((1, global_data.nN))
    d = {}
    for i in range(int(global_data.st)):
        if i == 0:
            temps = g.get_temps(True, temps)[0]
            g.calculate_matrix()
            g.to_file()
        else:
            g.calculate_matrix()
        g.matrix_aggregation()
        temps = g.solve_ode(temps)
        if i % 1000 == 0:
            print(f"poszlo {i}")
        if i % 5000 == 0:
            d[f'Time {i}'] = temps
            g.showMeshPlot()
    with open("Temperatures.csv", "w") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(d.keys())
        writer.writerows(zip(*d.values()))


def run():
    """
    Function calls calculate grid and runs simulation
    """
    g = FEM_Grid(True, 0)
    made(g)
