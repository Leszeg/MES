import matplotlib.pyplot as plt
import numpy as np
import xlwt
from numpy import zeros
from numpy.linalg import solve

from MES import Node as n, Element as e
from MES.Data import global_data


# Przechowuje listę elementów i listę węzłów potrzebne do stworzenia siatki
class FEM_Grid(object):
    """
    Class represents FEM grid

    Attributes
    ----------
    nodes : list
        List of nodes in gird

    elements : list
        List of elements in grid

    temperature_of_nodes : list
        List of nodes temperatures

    H_global : ndarray
        H_global matrix

    C_global : ndarray
        C_global matrix

    P_global : ndarray
        P_global matrix

    """

    def __init__(self, is_first: bool, temps):
        """
        Constructs all the necessary attributes for the FEM_grid object

        Parameters
        ----------
        is_first : bool
            Check if it is first create girid
        temps : list
            List of temperatures in nodes
        """

        d_x = global_data.B / (global_data.N_B - 1)
        d_y = global_data.H / (global_data.N_H - 1)

        self.nodes = []
        self.elements = []
        self.get_temps(is_first, temps)
        k = 0
        choice = True
        # Tworzenie współrzędnych węzłów
        for i1 in range(global_data.N_B):
            if i1 % 6 == 0 and i1 != 0:
                choice = not choice
            if not choice:
                for p in range(global_data.npm * 3):
                    pass
                    # Warunki odpowiadają za właściwe ustawienie warunku brzegowego(flaga BC) na krawędziach siatki
                    if p == 0:
                        self.nodes.append(n.Node(i1 * d_x, p * d_y, self.temperature_of_nodes[0][k], True))
                    elif p == 9 - 1:
                        self.nodes.append(n.Node(i1 * d_x, p * d_y, self.temperature_of_nodes[0][k], True))
                    else:
                        self.nodes.append(n.Node(i1 * d_x, p * d_y, self.temperature_of_nodes[0][k], False))
            else:
                for j1 in range(global_data.N_H):
                    # Warunki odpowiadają za właściwe ustawienie warunku brzegowego(flaga BC) na krawędziach siatki
                    if i1 == 0:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    elif j1 == 0:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    elif j1 == global_data.N_H - 1:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    elif i1 == global_data.N_B - 1:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    elif i1 / 5 == 1 and j1 > 7 or i1 / 17 == 1 and j1 > 7:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    elif i1 / 12 == 1 and j1 > 7 or i1 / 24 == 1 and j1 > 7:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    elif i1 / 29 == 1 and j1 > 7 or i1 / 36 == 1 and j1 > 7:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], True))
                    else:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], False))
            k += 1

        # Tworzenie elementów
        # Pętla for ma dodatek '+ global_data.N_B - 1' ponieważ przy tworzeniu siatki
        # elementy muszą być odpowiednio numerowane i będą dodatkowe 'puste przebiegi'
        # aby zachować odpowiednią numeracja przy końcach i początkach kolumn siatki
        tmp = [0, global_data.N_H, global_data.N_H + 1, 1, 0]
        choice = True
        for x2 in range(global_data.N_B + 2):
            if x2 / 5 == 1 and x2 != 1 or x2 / 12 == 1 or x2 / 17 == 1 or x2 / 24 == 1 or x2 / 29 == 1 or x2 / 36 == 1:
                choice = not choice
            if not choice:
                for x in range(7):
                    if x == 7 - 1:
                        tmp[0] += 1
                        tmp[1] += 1
                        tmp[2] += 1
                        tmp[3] += 1
                        tmp[4] += 1
                        break
                    ID = []
                    ID.append(tmp[4])
                    ID.append(ID[0] + global_data.N_H)
                    ID.append(ID[1] + 1)
                    ID.append(ID[0] + 1)
                    print(tmp)
                    print(ID)
                    nod = [self.nodes[tmp[0]], self.nodes[tmp[1]], self.nodes[tmp[2]], self.nodes[tmp[3]]]
                    self.elements.append(e.Element(ID, nod))
                    tmp[0] += 1
                    tmp[1] += 1
                    tmp[2] += 1
                    tmp[3] += 1
                    tmp[4] += 1
                print(f'koniec kolumny {x2}')
            else:
                for x1 in range(global_data.N_H):
                    if x1 == global_data.N_H - 1:
                        tmp[0] += 1
                        tmp[1] += 1
                        tmp[2] += 1
                        tmp[3] += 1
                        tmp[4] += 1
                        break
                    ID = []
                    ID.append(tmp[4])
                    ID.append(ID[0] + global_data.N_H)
                    ID.append(ID[1] + 1)
                    ID.append(ID[0] + 1)
                    print(tmp)
                    print(ID)
                    nod = [self.nodes[tmp[0]], self.nodes[tmp[1]], self.nodes[tmp[2]], self.nodes[tmp[3]]]
                    self.elements.append(e.Element(ID, nod))
                    tmp[0] += 1
                    tmp[1] += 1
                    tmp[2] += 1
                    tmp[3] += 1
                    tmp[4] += 1
                print(f'koniec kolumny {x2}')
        self.H_global, self.C_global, self.P_global = self.matrix_aggregation()

    def calculate_matrix(self):
        for j in range(len(self.elements)):
            self.elements[j].H_matrix()
            self.elements[j].C_matrix()
            self.elements[j].boundary_condition()

    def solve_ode(self, temps):
        Hz = self.H_global + (self.C_global / global_data.sst)
        x = -np.dot(self.C_global / global_data.sst, temps.T).T
        Pz = -(self.P_global + x)
        x = solve(Hz, Pz.T).T
        self.get_temps(False, x)
        return x

    def plot_grid(self):
        """
        Function plot grid elements and nodes.
        Blue - simple node
        Red - boundary condition
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for element in self.elements:
            for i in range(4):
                if element.nodes[i].bc == 0:
                    plt.scatter(element.x[i], element.y[i], color='blue')
                    if len(self.elements) < 10:
                        plt.annotate(element.nodes_ID[i], (element.x[i], element.y[i]))
                if element.nodes[i].bc == 1:
                    plt.scatter(element.x[i], element.y[i], color='red')
                    if len(self.elements) < 10:
                        plt.annotate(element.nodes_ID[i], (element.x[i], element.y[i]))

        ax.grid(which='both')
        plt.show()
        # plt.savefig("plot.png")

    def showMeshPlot(self):
        x = []
        y = []
        for node in self.nodes:
            x.append(node.x)
            y.append(node.y)
        plt.plot(x, y, ".k")
        plt.show()

    def matrix_aggregation(self):
        """
        Function aggregate local matrices to global
        Parameters
        ----------
        H_locals : ndarray
        C_locals : ndarray
        P_locals : ndarray

        Returns
        -------
        tuple
            Returns global matrices
        """
        r = global_data.N_B * global_data.N_H
        Pg = zeros(r, float)
        Hg = zeros((r, r), float)
        Cg = zeros((r, r), float)
        for i in range(len(self.elements)):
            h1 = self.elements[i].H_matrix_for_element
            c1 = self.elements[i].C_matrix_for_element
            for j in range(4):
                p1 = self.elements[i].P_matrix_for_element
                Pg[self.elements[i].nodes_ID[j]] += p1[j]
                for k in range(4):
                    Hg[self.elements[i].nodes_ID[j]][self.elements[i].nodes_ID[k]] += h1[j][k]
                    Cg[self.elements[i].nodes_ID[j]][self.elements[i].nodes_ID[k]] += c1[j][k]
        return Hg, Cg, Pg

    def to_file(self):
        book = xlwt.Workbook(encoding="utf-8")
        sheet = []
        counter = 0
        for omg in range(global_data.nE):
            sheet.append(book.add_sheet(f"ELement_{omg}"))
            element = self.elements[omg]

            sheet[omg].write(0, 0, "ID węzłów")
            for i in range(4):
                sheet[omg].write(1, i, element.nodes_ID[i])

            sheet[omg].write(3, 0, "Współrzędne X")
            for i in range(4):
                sheet[omg].write(4, i, element.x[i])

            sheet[omg].write(6, 0, "Współrzędne Y")
            for i in range(4):
                sheet[omg].write(7, i, element.y[i])

            sheet[omg].write(9, 0, "Ilość punktów całkowania")
            sheet[omg].write(10, 0, global_data.ip)

            sheet[omg].write(12, 0, "Współrzędne ksi punktów całkowania")
            for i in range(global_data.ip):
                sheet[omg].write(13, i, element.ksi[i])

            sheet[omg].write(15, 0, "Współrzędne eta punktów całkowania")
            for i in range(global_data.ip):
                sheet[omg].write(16, i, element.eta[i])

            sheet[omg].write(18, 0, "Macierz Jacobiego")
            for j in range(global_data.ip):
                counter = 19 + j
                for i in range(4):
                    sheet[omg].write(counter, i, element.jacob_matrix[j][i])

            sheet[omg].write(counter + 2, 0, "dN_dKsi")
            counter += 2
            for j in range(4):
                counter += 1
                for i in range(global_data.ip):
                    sheet[omg].write(counter, i, element.dN_dKsi[j][i])

            sheet[omg].write(counter + 2, 0, "dN_dEta")
            counter += 2
            for j in range(4):
                counter += 1
                for i in range(global_data.ip):
                    sheet[omg].write(counter, i, element.dN_dEta[j][i])

            sheet[omg].write(counter + 2, 0, "Wyznacznik macierzy Jacobiego")
            counter += 3
            for i in range(global_data.ip):
                sheet[omg].write(counter, i, element.determinant[i])

            sheet[omg].write(counter + 2, 0, "Odwrócona macierz Jacobiego")
            counter += 2
            for j in range(global_data.ip):
                counter += 1
                for i in range(4):
                    sheet[omg].write(counter, i, element.inv_jac[j][i])
            tmp2 = counter
            for i in range(global_data.ip):
                sheet[omg].write(counter + 2, 0, f"Macierz H dla {i} punktu całkowania")
                counter += 2
                H = element.H_matrix_for_ip[i]
                for j in range(4):
                    counter += 1
                    for k in range(4):
                        sheet[omg].write(counter, k, H[j][k])

            for i in range(global_data.ip):
                sheet[omg].write(tmp2 + 2, 6, f"Macierz C dla {i} punktu całkowania")
                tmp2 += 2
                C = element.C_matrix_for_ip[i]
                for j in range(4):
                    tmp2 += 1
                    for k in range(4):
                        sheet[omg].write(tmp2, 6 + k, C[j][k])

            sheet[omg].write(counter + 2, 0, "Macierz H dla elementu")
            counter += 2
            for i in range(4):
                counter += 1
                for j in range(4):
                    sheet[omg].write(counter, j, element.H_matrix_for_element[i][j])

            sheet[omg].write(tmp2 + 2, 6, "Macierz C dla elementu")
            tmp2 += 2
            for i in range(4):
                tmp2 += 1
                for j in range(4):
                    sheet[omg].write(tmp2, j + 6, element.C_matrix_for_element[i][j])
            xD = omg

        book.save("MES.xls")

    def get_temps(self, is_first, temps):
        if is_first:
            self.temperature_of_nodes = np.full((1, global_data.N_B * global_data.N_H), global_data.it)
            return self.temperature_of_nodes
        else:
            self.temperature_of_nodes = temps
            return self.temperature_of_nodes
