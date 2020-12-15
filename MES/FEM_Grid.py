import xlwt
from MES.Data import global_data
from MES import Node as n, Element as e
from numpy import arange, zeros, round
import matplotlib.pyplot as plt
import numpy as np


# Przechowuje listę elementów i listę węzłów potrzebne do stworzenia siatki
class FEM_Grid:

    def __init__(self, is_first: bool, temps):
        # Odległości między węzłami na odpowiednich osiach
        d_x = global_data.B / (global_data.N_B - 1)
        d_y = global_data.H / (global_data.N_H - 1)
        nE = global_data.nE
        # Inicjowanie pustych tablic
        self.nodes = []
        self.elements = []
        if is_first:
            for i in range(global_data.nN):
                self.temperature_of_nodes = np.full((1, global_data.N_B * global_data.N_H), global_data.it)
        else:
            self.temperature_of_nodes = temps
        k = 0
        # Tworzenie współrzędnych węzłów
        for i1 in range(global_data.N_B):
            for j1 in range(global_data.N_H):
                # Warunki odpowiadają za właściwe ustawienie warunku brzegowego(flaga BC) na krawędziach siatki
                if i1 == 0:
                    self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                elif j1 == 0:
                    self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                elif j1 == global_data.N_H - 1:
                    self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                elif i1 == global_data.N_B - 1:
                    self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 1))
                else:
                    self.nodes.append(n.Node(i1 * d_x, j1 * d_y, self.temperature_of_nodes[0][k], 0))
                k += 1

        # Tworzenie elementów
        # Pętla for ma dodatek '+ int((global_data.nW) / 2)' ponieważ przy tworzeniu siatki
        # elementy muszą być odpowiednio numerowane i będą dodatkowe 'puste przebiegi'
        # aby zachować odpowiednią numeracja przy końcach i początkach kolumn siatki
        tmp = [0, global_data.N_H, global_data.N_H + 1, 1, 0]
        column_end = 0
        for i in range(global_data.nE + global_data.N_B - 1):
            if column_end < global_data.N_H - 1:
                ID = []
                ID.append(tmp[4])
                ID.append(ID[0] + global_data.N_H)
                ID.append(ID[1] + 1)
                ID.append(ID[0] + 1)
                nod = [self.nodes[tmp[0]], self.nodes[tmp[1]], self.nodes[tmp[2]], self.nodes[tmp[3]]]
                self.elements.append(e.Element(global_data.pc, ID, nod))
                column_end += 1
                tmp[0] += 1
                tmp[1] += 1
                tmp[2] += 1
                tmp[3] += 1
                tmp[4] += 1
            else:
                tmp[0] += 1
                tmp[1] += 1
                tmp[2] += 1
                tmp[3] += 1
                tmp[4] += 1
                column_end = 0

        H = []
        C = []
        P = []
        BCH = []
        for i in self.elements:
            H.append(i.H_matrix_for_element)
            C.append(i.C_matrix_for_element)
            P.append(i.P_matrix_for_element)
            BCH.append(i.BCH_matrix_for_element)

        self.H_global = self.matrix_global(H)
        self.C_global = self.matrix_global(C)
        self.P_global = self.P_global(P)
        self.BCH_global = self.matrix_global(BCH)

    def plot_grid(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for element in self.elements:
            for i in range(4):
                if element.nodes[i].bc == 0:
                    plt.scatter(element.x[i], element.y[i], color='blue')
                    if (len(self.elements) < 10):
                        plt.annotate(element.nodes_ID[i], (element.x[i], element.y[i]))
                if element.nodes[i].bc == 1:
                    plt.scatter(element.x[i], element.y[i], color='red')
                    if (len(self.elements) < 10):
                        plt.annotate(element.nodes_ID[i], (element.x[i], element.y[i]))

        ax.grid(which='both')
        plt.show()

    def matrix_global(self, H_locals):
        r = global_data.N_B * global_data.N_H
        Hg = zeros((r, r), float)
        for i in range(len(self.elements)):
            h1 = H_locals[i]
            for j in range(4):
                for k in range(4):
                    Hg[self.elements[i].nodes_ID[j]][self.elements[i].nodes_ID[k]] += h1[j][k]
        return Hg

    def P_global(self, P_locals):
        r = global_data.N_B * global_data.N_H
        Pg = zeros((r), float)
        for i in range(len(self.elements)):
            for j in range(4):
                Pg[self.elements[i].nodes_ID[j]] += P_locals[i][j]
        return Pg

    def print_element_data(self, e):
        element = self.elements[e]
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

    def print_global_matrix(self):
        print("Globalna macierz H")
        print(self.H_global)

        print("Globlna macierz C")
        print(self.C_global)

    def to_file(self):
        book = xlwt.Workbook(encoding="utf-8")
        sheet = []

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
            sheet[omg].write(10, 0, element.integration_points_count)

            sheet[omg].write(12, 0, "Współrzędne ksi punktów całkowania")
            for i in range(global_data.pc):
                sheet[omg].write(13, i, element.ksi[i])

            sheet[omg].write(15, 0, "Współrzędne eta punktów całkowania")
            for i in range(global_data.pc):
                sheet[omg].write(16, i, element.eta[i])

            sheet[omg].write(18, 0, "Macierz jacobiego")
            for j in range(global_data.pc):
                tmp = 19 + j
                for i in range(4):
                    sheet[omg].write(tmp, i, element.jacob_matrix[j][i])

            sheet[omg].write(tmp + 2, 0, "dN_dKsi")
            tmp += 2
            for j in range(4):
                tmp += 1
                for i in range(global_data.pc):
                    sheet[omg].write(tmp, i, element.dN_dKsi[j][i])

            sheet[omg].write(tmp + 2, 0, "dN_dEta")
            tmp += 2
            for j in range(4):
                tmp += 1
                for i in range(global_data.pc):
                    sheet[omg].write(tmp, i, element.dN_dEta[j][i])

            sheet[omg].write(tmp + 2, 0, "Wyznacznik macierzy Jacobiego")
            tmp += 3
            for i in range(global_data.pc):
                sheet[omg].write(tmp, i, element.determinant[i])

            sheet[omg].write(tmp + 2, 0, "Odwrócona macierz Jacobiego")
            tmp += 2
            for j in range(global_data.pc):
                tmp += 1
                for i in range(4):
                    sheet[omg].write(tmp, i, element.inv_jac[j][i])
            tmp2 = tmp
            for i in range(global_data.pc):
                sheet[omg].write(tmp + 2, 0, f"Macierz H dla {i} punktu całkowania")
                tmp += 2
                H = element.H_matrix_for_ip[i]
                for j in range(4):
                    tmp += 1
                    for k in range(4):
                        sheet[omg].write(tmp, k, H[j][k])

            for i in range(global_data.pc):
                sheet[omg].write(tmp2 + 2, 6, f"Macierz C dla {i} punktu całkowania")
                tmp2 += 2
                C = element.C_matrix_for_ip[i]
                for j in range(4):
                    tmp2 += 1
                    for k in range(4):
                        sheet[omg].write(tmp2, 6 + k, C[j][k])

            sheet[omg].write(tmp + 2, 0, "Macierz H dla elementu")
            tmp += 2
            for i in range(4):
                tmp += 1
                for j in range(4):
                    sheet[omg].write(tmp, j, element.H_matrix_for_element[i][j])

            sheet[omg].write(tmp2 + 2, 6, "Macierz C dla elementu")
            tmp2 += 2
            for i in range(4):
                tmp2 += 1
                for j in range(4):
                    sheet[omg].write(tmp2, j + 6, element.C_matrix_for_element[i][j])
            xD = omg

        sheet.append(book.add_sheet("Macierz H globalna"))

        for j in range(global_data.N_B * global_data.N_H):
            for i in range(global_data.N_B * global_data.N_H):
                sheet[xD + 1].write(j, i, self.H_global[j][i])

        sheet.append(book.add_sheet("Macierz C globalna"))
        for j in range(global_data.N_B * global_data.N_H):
            for i in range(global_data.N_B * global_data.N_H):
                sheet[xD + 2].write(j, i, self.C_global[j][i])

        sheet.append(book.add_sheet("Macierz BCH globalna"))
        for j in range(global_data.N_B * global_data.N_H):
            for i in range(global_data.N_B * global_data.N_H):
                sheet[xD + 3].write(j, i, self.BCH_global[j][i])

        sheet.append(book.add_sheet("Macierz P globalna"))
        for j in range(global_data.N_B * global_data.N_H):
            sheet[xD + 4].write(0, j, self.P_global[j])

        book.save("MES.xls")
