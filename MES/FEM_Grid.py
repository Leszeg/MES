import xlwt
from MES import global_data
from numpy import arange, zeros, round
import matplotlib.pyplot as plt


# Przechowuje listę elementów i listę węzłów potrzebne do stworzenia siatki
class FEM_Grid:

    def __init__(self, E):
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
        for j in range(global_data.nE):
            element = self.ELEM[j]
            for i in range(4):
                if element.nodes[i].bc == 0:
                    plt.scatter(element.x[i], element.y[i], color='blue')
                    plt.annotate(element.nodes_ID[i], (element.x[i], element.y[i]))
                if element.nodes[i].bc == 1:
                    plt.scatter(element.x[i], element.y[i], color='red')
                    plt.annotate(element.nodes_ID[i], (element.x[i], element.y[i]))

        major_ticks_x = arange(0, 0.1, 0.0333)
        major_ticks_y = arange(0, 0.21, 0.05)
        ax.set_xticks(major_ticks_x)
        ax.set_yticks(major_ticks_y)
        ax.grid(which='both')
        plt.show()

    def matrix_global(self, H_locals):
        r = global_data.nW * global_data.nH
        Hg = zeros((r, r), float)
        for i in range(len(self.ELEM)):
            h1 = H_locals[i]
            for j in range(4):

                for k in range(4):
                    Hg[self.ELEM[i].nodes_ID[j]][self.ELEM[i].nodes_ID[k]] += h1[j][k]
        return Hg

    def P_global(self, P_locals):
        r = global_data.nW * global_data.nH
        Pg = zeros((r), float)
        for i in range(len(self.ELEM)):
            for j in range(4):
                Pg[self.ELEM[i].nodes_ID[j]] += P_locals[i][j]
        return Pg

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

    def to_file(self):
        book = xlwt.Workbook(encoding="utf-8")
        sheet = []

        for omg in range(global_data.nE):
            sheet.append(book.add_sheet(f"ELement_{omg}"))
            element = self.ELEM[omg]

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

        for j in range(global_data.nW * global_data.nH):
            for i in range(global_data.nW * global_data.nH):
                sheet[xD + 1].write(j, i, self.H_global[j][i])

        sheet.append(book.add_sheet("Macierz C globalna"))
        for j in range(global_data.nW * global_data.nH):
            for i in range(global_data.nW * global_data.nH):
                sheet[xD + 2].write(j, i, self.C_global[j][i])

        book.save("MES.xls")
