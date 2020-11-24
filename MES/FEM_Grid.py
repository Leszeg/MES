from MES import global_data
from numpy import arange, zeros, round
import matplotlib.pyplot as plt


# Przechowuje listę elementów i listę węzłów potrzebne do stworzenia siatki
class FEM_Grid:

    def __init__(self, N, E):
        self.ND = list(N)
        self.ELEM = list(E)
        H = []
        C = []
        for i in self.ELEM:
            H.append(i.H_matrix_for_element)
            C.append(i.C_matrix_for_element)

        self.H_global = self.matrix_global(H)
        self.C_global = self.matrix_global(C)

    def plot_grid(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        for i in range(global_data.nN):
            plt.scatter(self.ND[i].x, self.ND[i].y, color='blue')
            plt.annotate(i, (self.ND[i].x, self.ND[i].y))

        major_ticks_x = arange(0, 0.1, 0.0333)
        major_ticks_y = arange(0, 0.21, 0.05)
        ax.set_xticks(major_ticks_x)
        ax.set_yticks(major_ticks_y)
        ax.grid(which='both')
        plt.show()

    def matrix_global(self, H_locals):
        r = global_data.nW * global_data.nH
        Hg = zeros((r, r), float)
        for k in range(len(self.ELEM)):
            h1 = H_locals[k]
            for i in range(4):
                for j in range(4):
                    Hg[self.ELEM[k].nodes_ID[i]][self.ELEM[k].nodes_ID[j]] += h1[i][j]
        return Hg

    def print_element_data(self, e):
        element = self.ELEM[e]
        print(f"Informacje o {e} elemencie siatki\n")
        print("ID węzłów")
        print(element.nodes_ID)
        print("\n")
        print("Współrzędne X")
        print(element.x)
        print("\n")
        print("Współrzędne Y")
        print(element.y)
        print("\n")
        print("Ilość punktów całkowania")
        print(element.integration_points_count)
        print("\n")
        print("Współrzędne ksi punktów całkowania")
        print(element.ksi)
        print("\n")
        print("Współrzędne eta punktów całkowania")
        print(element.eta)
        print("\n")
        print("Macierz jacobiego")
        print(element.jacob_matrix)
        print("\n")
        print("dN_dKsi")
        print(round(element.dN_dKsi, 4))
        print("\n")
        print("dN_dEta")
        print(round(element.dN_dEta, 4))
        print("\n")
        print("Wyznacznik macierzy Jacobiego")
        print(round(element.determinant, 8))
        print("\n")
        print("Odwrócona macierz Jacobiego")
        print(element.inv_jac)
        print("\n")
        for i in range(element.integration_points_count):
            print(f"Macierze H dla {i} punktu całkowania")
            print(element.H_matrix_for_ip[i])
            print("\n")
        print("Macierz H dla całego elementu ")
        print(element.H_matrix_for_element)
        print("\n")
        for i in range(element.integration_points_count):
            print(f"Macierze C dla {i} punktu całkowania")
            print(element.C_matrix_for_ip[i])
            print("\n")
        print("Macierz C dla całego elementu")
        print(element.C_matrix_for_element)
        print("\n")

    def print_grid_data(self):
        print("GLOBAL DATA")
        print(f"W = {global_data.W}")
        print(f"H = {global_data.H}")
        print(f"nH = {global_data.nH}")
        print(f"nW = {global_data.nW}")
        print(f"k = {global_data.k}")
        print(f"ro = {global_data.ro}")
        print(f"c = {global_data.c}")
        print(f"t0 = {global_data.t0}")

    def print_global_matrix(self):
        print("Globalna macierz H")
        print(self.H_global)

        print("Globlna macierz C")
        print(self.C_global)
