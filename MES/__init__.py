import csv
from multiprocessing import Process
from MES.Data import global_data
from MES.FEM_Grid import FEM_Grid
from numpy import zeros, dot
from numpy.linalg import solve


def get_temps(g: FEM_Grid):
    x = g.temperature_of_nodes.T
    return x.T


temps = zeros((global_data.N_B * global_data.N_H, 1))
g = FEM_Grid(True, temps)


def made():
    """
    Function run simulation

    """
    temps = zeros((global_data.N_B * global_data.N_H, 1))
    d = {}
    for i in range(0, int(global_data.st), int(global_data.sst)):
        if i == 0:
            global g
            temps = get_temps(g)
            g.to_file()
        else:
            g = FEM_Grid(False, temps)
            temps = get_temps(g)
        d[f'Node {i}'] = temps[0]
        Hz = g.H_global + (g.C_global / global_data.sst)
        x = -dot(g.C_global / global_data.sst, temps.T).T
        Pz = -(g.P_global + x)
        temps = solve(Hz, Pz.T).T
    keys = sorted(d.keys())
    with open("test.csv", "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(keys)
        writer.writerows(zip(*[d[key] for key in keys]))


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
