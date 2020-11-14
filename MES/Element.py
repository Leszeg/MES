from math import sqrt
from MES.Node import Node
import numpy as np


class Element:

    def __init__(self, c, Id, nodes_c):
        self.nodes_ID = list(Id)

        # Dodane na potrzeby całkowania informacja o ilości węzłów w elemencie
        self.nodes_count = c

        self.nodes = nodes_c

        if c == 4:
            # Na sztywno dodaje współrzędne ksi i eta do testowania
            node_v = [-1 / sqrt(3), 1 / sqrt(3)]
            self.nodes[0].ksi = node_v[0]
            self.nodes[0].eta = node_v[0]

            self.nodes[1].ksi = node_v[1]
            self.nodes[1].eta = node_v[0]

            self.nodes[2].ksi = node_v[1]
            self.nodes[2].eta = node_v[1]

            self.nodes[3].ksi = node_v[0]
            self.nodes[3].eta = node_v[1]

        elif c == 9:
            # Na sztywno dodaje współrzędne ksi i eta do testowania
            node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]

            self.nodes[0].ksi = node_v[0]
            self.nodes[0].eta = node_v[0]

            self.nodes[1].ksi = node_v[1]
            self.nodes[1].eta = node_v[0]

            self.nodes[2].ksi = node_v[2]
            self.nodes[2].eta = node_v[0]

            self.nodes[3].ksi = node_v[0]
            self.nodes[3].eta = node_v[1]

            self.nodes[4].ksi = node_v[1]
            self.nodes[4].eta = node_v[1]

            self.nodes[5].ksi = node_v[2]
            self.nodes[5].eta = node_v[1]

            self.nodes[6].ksi = node_v[0]
            self.nodes[6].eta = node_v[2]

            self.nodes[7].ksi = node_v[1]
            self.nodes[7].eta = node_v[2]

            self.nodes[8].ksi = node_v[2]
            self.nodes[8].eta = node_v[2]

    def integral(self):
        result = 0
        fpc = []
        tmp = 0

        if self.nodes_count == 4:
            # Dwa punkty 2D
            # print('1. -2yx^2 + 2xy + 4')  # na zajęciach wyszło 16 dla 2 punktów

            # Wagi
            Ak = [1, 1]

            for j in range(2):
                for i in range(2):
                    # Tutaj implementuje wzór z którego liczę całkę
                    fpc.append(
                        -2 * self.nodes[tmp].eta * self.nodes[tmp].ksi ** 2 + 2 * self.nodes[tmp].ksi *
                        self.nodes[
                            tmp].eta + 4)
                    result = result + fpc[tmp] * Ak[i] * Ak[j]
                    tmp += 1
            return result

        elif int(self.nodes_count) == 9:
            # Trzy punkty 2D
            # print('2. -5yx^2 + 2xy^2 + 10')  # na zajęciach wyszło 40 dla 3 punktów

            # Wagi
            Ak = [5 / 9, 8 / 9, 5 / 9]

            for j in range(3):
                for i in range(3):
                    # Tutaj implementuje wzór z którego liczę całkę
                    fpc.append(-5 * self.nodes[tmp].eta * self.nodes[tmp].ksi ** 2 + 2 * self.nodes[tmp].ksi *
                               self.nodes[tmp].eta ** 2 + 10)
                    result = result + fpc[tmp] * Ak[i] * Ak[j]
                    tmp += 1
            return result

    def jacobian(self):
        # Współrzędne globalne elementu
        x = [0, 0.025, 0.025, 0]
        y = [0, 0, 0.025, 0.025]
        # for i in range(4):
        #     x.append(self.nodes[i].x)
        #     y.append(self.nodes[i].y)

        # Współrzędne lokalne elementu numeracja jak z siatką FEM
        ksi = [self.nodes[0].ksi, self.nodes[1].ksi, self.nodes[2].ksi, self.nodes[3].ksi]
        eta = [self.nodes[0].eta, self.nodes[1].eta, self.nodes[2].eta, self.nodes[3].eta]

        # Obliczam pochodne funkcji kształtu po ksi i eta
        dN_dKsi = np.zeros((4, 4), float)
        dN_dEta = np.zeros((4, 4), float)
        for i in range(0, self.nodes_count):
            for j in range(0, self.nodes_count):
                if j == 0:
                    dN_dKsi[j][i] = -0.25 * (1 - eta[i])
                    dN_dEta[j][i] = -0.25 * (1 - ksi[i])
                elif j == 1:
                    dN_dKsi[j][i] = 0.25 * (1 - eta[i])
                    dN_dEta[j][i] = -0.25 * (1 + ksi[i])
                elif j == 2:
                    dN_dKsi[j][i] = 0.25 * (1 + eta[i])
                    dN_dEta[j][i] = 0.25 * (1 + ksi[i])
                elif j == 3:
                    dN_dKsi[j][i] = -0.25 * (1 + eta[i])
                    dN_dEta[j][i] = 0.25 * (1 - ksi[i])

        jacobian = np.zeros((4, 4), float)
        determinant = []
        inv_jacobian = np.zeros((4, 4), float)
        m = dN_dEta
        dN_dEta = dN_dKsi
        dN_dKsi = m
        # Cztery jakobiany 2x2 ułożone wierszami
        for i in range(0, self.nodes_count):
            for j in range(0, self.nodes_count):
                for k in range(0, self.nodes_count):
                    if j == 0:
                        jacobian[i][j] += dN_dEta[k][j] * x[k]
                    elif j == 1:
                        jacobian[i][j] += dN_dKsi[k][j] * x[k]
                    elif j == 2:
                        jacobian[i][j] += dN_dEta[k][j] * y[k]
                    elif j == 3:
                        jacobian[i][j] += dN_dKsi[k][j] * y[k]

        for i in range(self.nodes_count):
            determinant.append(jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2])

        k = 0
        for i in range(self.nodes_count):
            for j in range(3, -1, -1):
                if k == 0 or k == 3:
                    inv_jacobian[i][k] = round(jacobian[i][j] / determinant[i],3)
                elif k == 1 or k == 2:
                    inv_jacobian[i][k] = round(jacobian[i][j] / determinant[i],3)
                k += 1
            k = 0
        print("///////////X////////////")
        print(x)
        print("///////////Y////////////")
        print(y)

        print("//////////ETA///////////")
        for i in range(self.nodes_count):
            print(eta[i])

        print("//////////KSI///////////")
        for i in range(self.nodes_count):
            print(ksi[i])

        # Pochodne w macierzy ułożone wierszami - pierwszy wiersz - dN1_dEta itd.
        print('////////// POCHODNA FUNKCJI KSZTAŁTU ETA //////////')
        for i in range(self.nodes_count):
            print(dN_dKsi[i])
        # Pochodne w macierzy ułożone wierszami - pierwszy wiersz - dN1_dKsi itd.
        print('////////// POCHODNA FUNKCJI KSZTAŁTU KSI //////////')
        for i in range(self.nodes_count):
            print(dN_dEta[i])

        print('////////// JACOBIAN //////////')
        for i in range(self.nodes_count):
            print(jacobian[i])

        print('///////// DETERMINANT //////////')
        for i in range(self.nodes_count):
            print(determinant[i])

        print('////////// INVERSE JACOBIAN //////////')
        for i in range(self.nodes_count):
            print(inv_jacobian[i])

        return jacobian, dN_dKsi, dN_dEta, determinant, inv_jacobian

    def H_matrix(self):
        data = self.jacobian()
        jacobian = data[0]
        dN_dEta = data[1]
        dN_dKsi = data[2]
        determinant = data[3]
        inv_jacobian = data[4]
        dN_dY = np.zeros((4, 4), float)
        dN_dX = np.zeros((4, 4), float)
        for i in range(4):
            for k in range(4):
                dN_dX[i][k] = (dN_dKsi[i][k] * inv_jacobian[0][0] + dN_dEta[i][k] * (inv_jacobian[0][1]))
                dN_dY[i][k] = (dN_dKsi[i][k] * inv_jacobian[0][2] + dN_dEta[i][k] * (inv_jacobian[0][3]))

        print("dN_dX")
        print(dN_dX)

        print("dN_dY")
        print(dN_dY)
        dNdX1 = np.zeros((4, 4), float)
        dNdY1 = np.zeros((4, 4), float)
        dNdX2 = np.zeros((4, 4), float)
        dNdY2 = np.zeros((4, 4), float)
        dNdX3 = np.zeros((4, 4), float)
        dNdY3 = np.zeros((4, 4), float)
        dNdX4 = np.zeros((4, 4), float)
        dNdY4 = np.zeros((4, 4), float)

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX1[i][j] = dN_dX[i][0] * dN_dX[j][0]
                dNdY1[i][j] = dN_dY[i][0] * dN_dY[j][0]

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX2[i][j] = dN_dX[i][1] * dN_dX[j][1]
                dNdY2[i][j] = dN_dY[i][1] * dN_dY[j][1]

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX3[i][j] = dN_dX[i][2] * dN_dX[j][2]
                dNdY3[i][j] = dN_dY[i][2] * dN_dY[j][2]

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX4[i][j] = dN_dX[i][3] * dN_dX[j][3]
                dNdY4[i][j] = dN_dY[i][3] * dN_dY[j][3]
        print(dNdX1)
        for i in range(0, 4):
            for j in range(0, 4):
                dNdX1[i][j] += dNdY1[i][j]
                dNdX2[i][j] += dNdY2[i][j]
                dNdX3[i][j] += dNdY3[i][j]
                dNdX4[i][j] += dNdY4[i][j]

        H = np.zeros((4, 4), float)
        for i in range(0, 4):
            for j in range(0, 4):
                H[i][j] = 30 * (dNdX1[i][j] + dNdX2[i][j] + dNdX3[i][j] + dNdX4[i][j]) * determinant[i]

        print("\n\n")

        print("////////// H //////////")
        print(H)
        return H
