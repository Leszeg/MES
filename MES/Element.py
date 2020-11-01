from math import sqrt
from MES.Node import Node
import numpy as np


class Element:

    def __init__(self, c):
        self.nodes = []
        self.nodes_count = c
        if c == 4:
            node_v = [-1 / sqrt(3), 1 / sqrt(3)]
            for j in range(2):
                for i in range(2):
                    self.nodes.append(Node(node_v[i], node_v[j], 0, 0))
        elif c == 9:
            node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
            for j in range(3):
                for i in range(3):
                    self.nodes.append(Node(node_v[i], node_v[j], 0, 0))

    def kons(self, ID):
        self.nodes_ID = ID

    def integral(self):
        result = 0
        fpc = []
        tmp = 0

        if self.nodes_count == 4:
            # Dwa punkty 2D
            print('1. -2yx^2 + 2xy + 4')  # na zajęciach wyszło 16 dla 2 punktów
            Ak = [1, 1]
            for j in range(2):
                for i in range(2):
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
            Ak = [5 / 9, 8 / 9, 5 / 9]
            for j in range(3):
                for i in range(3):
                    fpc.append(-5 * self.nodes[tmp].eta * self.nodes[tmp].ksi ** 2 + 2 * self.nodes[tmp].ksi *
                               self.nodes[tmp].eta ** 2 + 10)
                    result = result + fpc[tmp] * Ak[i] * Ak[j]
                    tmp += 1
            return result

    def jacobian(self):
        x = [0, 4, 4, 0]
        y = [0, 0, 4, 4]
        derivativeEta = np.zeros((4, 4), float)
        derivativeKsi = np.zeros((4, 4), float)
        eta = [-1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)]
        ksi = [-1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3)]

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
                    if inv_jacobian[k][i] == -0.0:
                        inv_jacobian[k][i] = 0.0
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
