from MES import global_data
from numpy import arange, zeros
import matplotlib.pyplot as plt


# Przechowuje listę elementów i listę węzłów potrzebne do stworzenia siatki
class FEM_Grid:

    def __init__(self, N, E):
        self.ND = list(N)
        self.ELEM = list(E)

    def print_grid(self):
        for i in range(global_data.nN):
            print(self.ND[i].x, " ", self.ND[i].y)

        for i in range(global_data.nE):
            print(f"Element {i}")
            for j in range(4):
                print(self.ELEM[i].nodes_ID[j])

    def plot_grid(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        for i in range(nN):
            plt.scatter(nodes[i].x, nodes[i].y, color='blue')
            plt.annotate(i, (nodes[i].x, nodes[i].y))

        major_ticks_x = arange(0, 0.1, 0.0333)
        major_ticks_y = arange(0, 0.21, 0.05)
        ax.set_xticks(major_ticks_x)
        ax.set_yticks(major_ticks_y)
        ax.grid(which='both')
        plt.show()

    def matrix_global(self, H_locals, elements):
        r = global_data.nW * global_data.nH
        Hg = zeros((r, r), float)
        h1 = []
        for k in range(len(elements)):
            h1 = H_locals[k]
            for i in range(4):
                for j in range(4):
                    Hg[elements[k].nodes_ID[i]][elements[k].nodes_ID[j]] += h1[i][j]