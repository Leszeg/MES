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
            print('1. -2yx^2 + 2xy + 4')  # na zajęciach wyszło 16 dla 2 punktów

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
            print('2. -5yx^2 + 2xy^2 + 10')  # na zajęciach wyszło 40 dla 3 punktów

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
        x = [0, 4, 4, 0]
        y = [0, 0, 6, 6]
        derivativeEta = np.zeros((4, 4), float)
        derivativeKsi = np.zeros((4, 4), float)

        for i in range(4):
            print(self.nodes[i].eta)

        for i in range(4):
            print(self.nodes[i].ksi)

        eta = [self.nodes[0].eta, self.nodes[1].eta, self.nodes[3].eta, self.nodes[2].eta]
        ksi = [self.nodes[0].ksi, self.nodes[1].ksi, self.nodes[3].ksi, self.nodes[2].ksi]

        for i in range(0, self.nodes_count):
            for j in range(0, self.nodes_count):
                if j == 0:
                    derivativeEta[i][j] = -0.25 * (1 - eta[i])
                    derivativeKsi[i][j] = (-0.25) * (1 - ksi[i])
                elif j == 1:
                    derivativeEta[i][j] = 0.25 * (1 - eta[i])
                    derivativeKsi[i][j] = -0.25 * (1 + ksi[i])
                elif j == 2:
                    derivativeEta[i][j] = 0.25 * (1 + eta[i])
                    derivativeKsi[i][j] = 0.25 * (1 + ksi[i])
                elif j == 3:
                    derivativeEta[i][j] = -0.25 * (1 + eta[i])
                    derivativeKsi[i][j] = 0.25 * (1 - ksi[i])

        jacobian = np.zeros((4, 4), float)
        determinant = []
        inv_jacobian = np.zeros((4, 4), float)

        for i in range(0, self.nodes_count):
            for j in range(0, self.nodes_count):
                for k in range(0, self.nodes_count):
                    if j == 0:
                        jacobian[j][i] += derivativeEta[j][k] * x[k]
                    elif j == 1:
                        jacobian[j][i] += derivativeKsi[j][k] * x[k]
                    elif j == 2:
                        jacobian[j][i] += derivativeEta[j][k] * y[k]
                    elif j == 3:
                        jacobian[j][i] += derivativeKsi[j][k] * y[k]

        for i in range(self.nodes_count):
            determinant.append(jacobian[0][i] * jacobian[3][i] - jacobian[1][i] * jacobian[2][i])

        k = 0
        for i in range(self.nodes_count):
            for j in range(3, -1, -1):
                if k == 0 or k == 3:
                    inv_jacobian[k][i] = jacobian[j][i] / determinant[i]
                elif k == 1 or k == 2:
                    inv_jacobian[k][i] = -jacobian[j][i] / determinant[i]
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

        print('////////// POCHODNA FUNKCJI KSZTAŁTU ETA //////////')
        for i in range(self.nodes_count):
            print(derivativeEta[i])

        print('////////// POCHODNA FUNKCJI KSZTAŁTU KSI //////////')
        for i in range(self.nodes_count):
            print(derivativeKsi[i])

        print('////////// JACOBIAN //////////')
        for i in range(self.nodes_count):
            print(jacobian[i])

        print('///////// DETERMINANT //////////')
        for i in range(self.nodes_count):
            print(determinant[i])

        print('////////// INVERSE JACOBIAN //////////')
        for i in range(self.nodes_count):
            print(inv_jacobian[i])

        dN_dx = np.zeros((4, 4), float)
        for i in range(0, self.nodes_count):
            k = 0
            for j in range(2):
                for l in range(2):
                    dN_dx[i][k] = ((1 / determinant[0]) * (
                            derivativeKsi[i][k] * jacobian[1][k] + derivativeEta[i][k] * (jacobian[3][k])))
                    k += 1

        print("////////// Pochodna funkcji kształtu X //////////")
        for i in range(self.nodes_count):
            print(dN_dx[i])
        return jacobian, derivativeEta, derivativeKsi, determinant

    def H_matrix(self):

        data = self.jacobian()
        jacobian = data[0]
        derivativeEta = data[1]
        derivativeKsi = data[2]
        determinant = data[3]

        jacSize = 2
        One = np.zeros((jacSize, jacSize), float)
        Two = np.zeros((jacSize, jacSize), float)
        Three = np.zeros((jacSize, jacSize), float)
        Four = np.zeros((jacSize, jacSize), float)
        k = 0

        for i in range(jacSize):
            for j in range(jacSize):
                One[j][i] = jacobian[k][0]
                Two[j][i] = jacobian[k][1]
                Three[j][i] = jacobian[k][2]
                Four[j][i] = jacobian[k][3]
                k += 1

        dNdX = np.zeros((4, 4), float)
        dNdY = np.zeros((4, 4), float)

        dNdX[0][0] = (1 / determinant[0]) * (One[1][0] * derivativeEta[0][0] + One[1][1] * derivativeEta[0][0])
        dNdX[0][1] = (1 / determinant[0]) * (One[1][0] * derivativeEta[0][1] + One[1][1] * derivativeEta[0][1])
        dNdX[0][2] = (1 / determinant[0]) * (One[1][0] * derivativeEta[0][2] + One[1][1] * derivativeEta[0][2])
        dNdX[0][3] = (1 / determinant[0]) * (One[1][0] * derivativeEta[0][3] + One[1][1] * derivativeEta[0][3])

        dNdX[1][0] = (1 / determinant[1]) * (Two[1][0] * derivativeEta[1][0] + Two[1][1] * derivativeEta[1][0])
        dNdX[1][1] = (1 / determinant[1]) * (Two[1][0] * derivativeEta[1][1] + Two[1][1] * derivativeEta[1][1])
        dNdX[1][2] = (1 / determinant[1]) * (Two[1][0] * derivativeEta[1][2] + Two[1][1] * derivativeEta[1][2])
        dNdX[1][3] = (1 / determinant[1]) * (Two[1][0] * derivativeEta[1][3] + Two[1][1] * derivativeEta[1][3])

        dNdX[2][0] = (1 / determinant[2]) * (
                Three[1][0] * derivativeEta[2][0] + Three[1][1] * derivativeEta[2][0])
        dNdX[2][1] = (1 / determinant[2]) * (
                Three[1][0] * derivativeEta[2][1] + Three[1][1] * derivativeEta[2][1])
        dNdX[2][2] = (1 / determinant[2]) * (
                Three[1][0] * derivativeEta[2][2] + Three[1][1] * derivativeEta[2][2])
        dNdX[2][3] = (1 / determinant[2]) * (
                Three[1][0] * derivativeEta[2][3] + Three[1][1] * derivativeEta[2][3])

        dNdX[3][0] = (1 / determinant[3]) * (Four[1][0] * derivativeEta[3][0] + Four[1][1] * derivativeEta[3][0])
        dNdX[3][1] = (1 / determinant[3]) * (Four[1][0] * derivativeEta[3][1] + Four[1][1] * derivativeEta[3][1])
        dNdX[3][2] = (1 / determinant[3]) * (Four[1][0] * derivativeEta[3][2] + Four[1][1] * derivativeEta[3][2])
        dNdX[3][3] = (1 / determinant[3]) * (Four[1][0] * derivativeEta[3][3] + Four[1][1] * derivativeEta[3][3])
        ##################################################################################################
        dNdY[0][0] = (1 / determinant[0]) * (One[0][0] * derivativeKsi[0][0] + One[0][1] * derivativeKsi[0][0])
        dNdY[0][1] = (1 / determinant[0]) * (One[0][0] * derivativeKsi[0][1] + One[0][1] * derivativeKsi[0][1])
        dNdY[0][2] = (1 / determinant[0]) * (One[0][0] * derivativeKsi[0][2] + One[0][1] * derivativeKsi[0][2])
        dNdY[0][3] = (1 / determinant[0]) * (One[0][0] * derivativeKsi[0][3] + One[0][1] * derivativeKsi[0][3])

        dNdY[1][0] = (1 / determinant[1]) * (Two[0][0] * derivativeKsi[1][0] + Two[0][1] * derivativeKsi[1][0])
        dNdY[1][1] = (1 / determinant[1]) * (Two[0][0] * derivativeKsi[1][1] + Two[0][1] * derivativeKsi[1][1])
        dNdY[1][2] = (1 / determinant[1]) * (Two[0][0] * derivativeKsi[1][2] + Two[0][1] * derivativeKsi[1][2])
        dNdY[1][3] = (1 / determinant[1]) * (Two[0][0] * derivativeKsi[1][3] + Two[0][1] * derivativeKsi[1][3])

        dNdY[2][0] = (1 / determinant[2]) * (
                Three[0][0] * derivativeKsi[2][0] + Three[0][1] * derivativeKsi[2][0])
        dNdY[2][1] = (1 / determinant[2]) * (
                Three[0][0] * derivativeKsi[2][1] + Three[0][1] * derivativeKsi[2][1])
        dNdY[2][2] = (1 / determinant[2]) * (
                Three[0][0] * derivativeKsi[2][2] + Three[0][1] * derivativeKsi[2][2])
        dNdY[2][3] = (1 / determinant[2]) * (
                Three[0][0] * derivativeKsi[2][3] + Three[0][1] * derivativeKsi[2][3])

        dNdY[3][0] = (1 / determinant[3]) * (Four[0][0] * derivativeKsi[3][0] + Four[0][1] * derivativeKsi[3][0])
        dNdY[3][1] = (1 / determinant[3]) * (Four[0][0] * derivativeKsi[3][1] + Four[0][1] * derivativeKsi[3][1])
        dNdY[3][2] = (1 / determinant[3]) * (Four[0][0] * derivativeKsi[3][2] + Four[0][1] * derivativeKsi[3][2])
        dNdY[3][3] = (1 / determinant[3]) * (Four[0][0] * derivativeKsi[3][3] + Four[0][1] * derivativeKsi[3][3])

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
                dNdX1[i][j] = dNdX[0][i] * dNdX[0][j]
                dNdY1[i][j] = dNdY[0][i] * dNdY[0][j]

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX2[i][j] = dNdX[1][i] * dNdX[1][j]
                dNdY2[i][j] = dNdY[1][i] * dNdY[1][j]

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX3[i][j] = dNdX[2][i] * dNdX[2][j]
                dNdY3[i][j] = dNdY[2][i] * dNdY[2][j]

        for i in range(0, 4):
            for j in range(0, 4):
                dNdX4[i][j] = dNdX[3][i] * dNdX[3][j]
                dNdY4[i][j] = dNdY[3][i] * dNdY[3][j]

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
        print(H[0])
        print(H[1])
        print(H[2])
        print(H[3])