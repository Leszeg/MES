import matplotlib.pyplot as plt
import matplotlib.collections
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
        self.temperature_of_nodes = self.temperature_of_nodes = np.full((1, global_data.nN), global_data.it)

        k = 0
        choice = True
        # Tworzenie współrzędnych węzłów
        for i1 in range(global_data.N_B):
            if i1 % 6 == 0 and i1 != 0:
                choice = not choice
            if not choice:
                for p in range(global_data.npm * 3):
                    # Warunki odpowiadają za właściwe ustawienie warunku brzegowego(flaga BC) na krawędziach siatki
                    if p == 0:
                        self.temperature_of_nodes[0][k] = 47
                        self.nodes.append(n.Node(i1 * d_x, p * d_y, self.temperature_of_nodes[0][k], 0))  # 47
                    elif p == 9 - 1:
                        self.nodes.append(n.Node(i1 * d_x, p * d_y, self.temperature_of_nodes[0][k], 1))
                    else:
                        self.nodes.append(n.Node(i1 * d_x, p * d_y, self.temperature_of_nodes[0][k], 0))
                    k += 1
            else:
                for j1 in range(global_data.N_H):
                    # Warunki odpowiadają za właściwe ustawienie warunku brzegowego(flaga BC) na krawędziach siatki
                    if i1 == 0 and j1 == 0:
                        self.temperature_of_nodes[0][k] = 47
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))  # 47
                    elif i1 == 0:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    elif i1 == global_data.N_B - 1:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    elif j1 == 0:
                        self.temperature_of_nodes[0][k] = 47
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 0))  # 47
                    elif j1 == global_data.N_H - 1:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    elif i1 == global_data.N_B - 1:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    elif i1 / 5 == 1 and j1 > 7 or i1 / 17 == 1 and j1 > 7:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    elif i1 / 12 == 1 and j1 > 7 or i1 / 24 == 1 and j1 > 7:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    elif i1 / 29 == 1 and j1 > 7 or i1 / 36 == 1 and j1 > 7:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                    else:
                        self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 0))
                    k += 1

        self.get_temps(is_first, temps)
        # Tworzenie elementów
        # Pętla for ma dodatek '+ global_data.N_B - 1' ponieważ przy tworzeniu siatki
        # elementy muszą być odpowiednio numerowane i będą dodatkowe 'puste przebiegi'
        # aby zachować odpowiednią numeracja przy końcach i początkach kolumn siatki
        tmp = [0, global_data.N_H, global_data.N_H + 1, 1, 0]
        choice = True
        k = 0
        for x2 in range(global_data.N_B - 1):
            if x2 / 5 == 1 or x2 / 12 == 1 or x2 / 17 == 1 or x2 / 24 == 1 or x2 / 29 == 1 or x2 / 36 == 1:
                choice = not choice
            if not choice:
                for o in range(9):
                    if o == 8:
                        if k == 0 or k == 7 or k == 14:
                            tmp[0] += 13
                            tmp[1] += 1
                            tmp[2] += 1
                            tmp[3] += 13
                            tmp[4] += 1
                        elif k == 6 or k == 13 or k == 20:
                            tmp[0] += 1
                            tmp[1] += 13
                            tmp[2] += 13
                            tmp[3] += 1
                            tmp[4] += 1

                        else:
                            tmp[0] += 1
                            tmp[1] += 1
                            tmp[2] += 1
                            tmp[3] += 1
                            tmp[4] += 1
                        k = k + 1
                        break
                    ID = []
                    ID.append(tmp[0])
                    ID.append(tmp[1])
                    ID.append(tmp[2])
                    ID.append(tmp[3])
                    nod = [self.nodes[tmp[0]], self.nodes[tmp[1]], self.nodes[tmp[2]], self.nodes[tmp[3]]]

                    if x2 == 0:
                        self.elements.append(e.Element(ID, nod, False))
                    elif x2 < 3 and x2 != 0:
                        self.elements.append(e.Element(ID, nod, True))
                    else:
                        self.elements.append(e.Element(ID, nod, False))
                    tmp[0] += 1
                    tmp[1] += 1
                    tmp[2] += 1
                    tmp[3] += 1
                    tmp[4] += 1
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
                    ID.append(tmp[0])
                    ID.append(tmp[1])
                    ID.append(tmp[2])
                    ID.append(tmp[3])
                    nod = [self.nodes[tmp[0]], self.nodes[tmp[1]], self.nodes[tmp[2]], self.nodes[tmp[3]]]
                    if x1 == 0:
                        self.elements.append(e.Element(ID, nod, False))
                    elif x1 < 3 and x1 != 0:
                        self.elements.append(e.Element(ID, nod, True))
                    else:
                        self.elements.append(e.Element(ID, nod, False))
                    tmp[0] += 1
                    tmp[1] += 1
                    tmp[2] += 1
                    tmp[3] += 1
                    tmp[4] += 1
        self.H_global, self.C_global, self.P_global = self.matrix_aggregation()

    def calculate_matrix(self):
        for j in range(len(self.elements)):
            self.elements[j].H_matrix()
            self.elements[j].C_matrix()
            self.elements[j].boundary_condition()

    def solve_ode(self, temps):

        Hz = self.H_global + (self.C_global / 0.002)
        x = -np.dot(self.C_global / 0.002, temps.T).T
        Pz = -(self.P_global + x)
        x = solve(Hz, Pz.T).T
        for i in range(global_data.nN):
            if temps[i] == 47:
                x[i] = 47
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
        li = []
        li2 = []
        for node in self.nodes:
            x.append(node.x)
            y.append(node.y)
        for i in range(global_data.nE):
            li.append(self.elements[i].nodes_ID)
            sr = self.temperature_of_nodes[li[i][0]] + self.temperature_of_nodes[li[i][1]] + self.temperature_of_nodes[
                li[i][2]] + \
                 self.temperature_of_nodes[li[i][3]]
            sr = sr / 4
            li2.append(sr)
        values = np.array(li2)
        elem = np.array(li)

        def quatplot(x, y, quatrangles, values, ax=None, **kwargs):
            if not ax: ax = plt.gca()
            xy = np.c_[x, y]
            verts = xy[quatrangles]
            pc = matplotlib.collections.PolyCollection(verts, **kwargs)
            pc.set_array(values)
            ax.add_collection(pc)
            ax.autoscale()
            return pc

        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        self.w = quatplot(x, y, np.asarray(elem), values, ax=ax,
                          edgecolor="crimson", cmap="rainbow")
        fig.colorbar(w, ax=ax)
        ax.plot(x, y, marker=",", ls="", color="crimson")

        ax.set(title='This is the plot for: quad', xlabel='X Axis', ylabel='Y Axis')

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
        r = global_data.nN
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
            return self.temperature_of_nodes
        else:
            self.temperature_of_nodes = temps
            return self.temperature_of_nodes
