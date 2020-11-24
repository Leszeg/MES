from math import sqrt
import numpy as np
from MES.__init__ import global_data
from MES import IntegrationPoint as IP


class Element:

    def __init__(self, c, Id, nodes_c):
        # Lista z odpowiednio ułożonymi ID węzłów do wydrukowania siatki
        self.nodes_ID = list(Id)

        # Ilość punktów całkowania w elemencie
        self.integration_points_count = int(c)

        # Współrzędne globalne elementu
        self.x = []
        self.y = []
        for i in range(4):
            self.x.append(nodes_c[i].x)
            self.y.append(nodes_c[i].y)

        # Deklaracja tablicy punktów całkowania
        integration_points = []

        # Odpowiednie współrzędne lokalne w zależności od ilości punktów
        if self.integration_points_count == 1:
            node_v = [1, 2, 3]
        if self.integration_points_count == 4:
            node_v = [-1 / sqrt(3), 1 / sqrt(3)]

        elif self.integration_points_count == 9:
            node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]

        elif self.integration_points_count == 16:
            node_v = [-0.861136, -0.339981, 0.339981, 0.861136]

        # Tworzymy tablicę punktów całkowania
        k = 0
        t = 0
        for i in range(self.integration_points_count):
            integration_points.append(IP.IntegrationPoint())
            integration_points[i].ksi = node_v[k]
            integration_points[i].eta = node_v[t]
            k += 1
            if k % len(node_v) == 0:
                t += 1
                k = 0

        self.ksi = []
        self.eta = []
        for i in range(self.integration_points_count):
            self.ksi.append(integration_points[i].ksi)
            self.eta.append(integration_points[i].eta)

        data = self.jacobian()
        self.jacob_matrix = data[0]
        self.dN_dEta = data[1]
        self.dN_dKsi = data[2]
        self.determinant = data[3]
        self.inv_jac = data[4]

        data = self.H_matrix()
        # Lista lokalnych macierzy H dla każdego punktu całkowania
        self.H_matrix_for_ip = data[0]

        # Macierz H dla elementu
        self.H_matrix_for_element = data[1]

        data = self.C_matrix()
        # Lista lokalnych macierzy H dla każdego punktu całkowania
        self.C_matrix_for_ip = data[0]

        # Macierz H dla elementu
        self.C_matrix_for_element = data[1]

    def integral(self, H_matrix):
        result = []
        k = 0

        # Odpowiednie wagi w zależności od ilości punktów całkowania
        if self.integration_points_count == 4:
            # Wagi
            Ak = [1, 1]

        elif int(self.integration_points_count) == 9:
            # Wagi
            Ak = [5 / 9, 8 / 9, 5 / 9]

        elif int(self.integration_points_count) == 16:
            # Wagi
            Ak = [0.347855, 0.652145, 0.652145, 0.347855]

        # Mnożymy przez wagi
        for j in range(len(Ak)):
            for i in range(len(Ak)):
                result.append(H_matrix[k] * Ak[i] * Ak[j])
                k += 1
        return result

    def jacobian(self):
        # Obliczam pochodne funkcji kształtu po ksi i eta
        dN_dKsi = np.zeros((4, self.integration_points_count), float)
        dN_dEta = np.zeros((4, self.integration_points_count), float)
        for i in range(0, self.integration_points_count):
            for j in range(0, self.integration_points_count):
                if j == 0:
                    dN_dKsi[j][i] = -0.25 * (1 - self.eta[i])
                    dN_dEta[j][i] = -0.25 * (1 - self.ksi[i])
                elif j == 1:
                    dN_dKsi[j][i] = 0.25 * (1 - self.eta[i])
                    dN_dEta[j][i] = -0.25 * (1 + self.ksi[i])
                elif j == 2:
                    dN_dKsi[j][i] = 0.25 * (1 + self.eta[i])
                    dN_dEta[j][i] = 0.25 * (1 + self.ksi[i])
                elif j == 3:
                    dN_dKsi[j][i] = -0.25 * (1 + self.eta[i])
                    dN_dEta[j][i] = 0.25 * (1 - self.ksi[i])

        jacobian = np.zeros((self.integration_points_count, 4), float)
        determinant = []
        inv_jacobian = np.zeros((self.integration_points_count, 4), float)

        # Cztery jakobiany 2x2 ułożone wierszami
        sum = 0
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for i in range(0, self.integration_points_count):
            for j in range(0, 4):
                sum += dN_dKsi[j][i] * self.x[j]
                sum1 += dN_dKsi[j][i] * self.y[j]
                sum2 += dN_dEta[j][i] * self.x[j]
                sum3 += dN_dEta[j][i] * self.y[j]
            jacobian[i][0] = sum
            jacobian[i][1] = sum1
            jacobian[i][2] = sum2
            jacobian[i][3] = sum3
            sum = 0
            sum1 = 0
            sum2 = 0
            sum3 = 0

        for i in range(self.integration_points_count):
            determinant.append(jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2])

        # zamiana wiersza jakobiany na macierz 2x2
        tmp_jacobian = []
        for i in range(self.integration_points_count):
            tmp_jacobian.append(np.zeros((2, 2), float))

        # Obliczanie odwrotnego jakobianu
        t = 0
        for k in range(self.integration_points_count):
            tmp = tmp_jacobian[k]
            for i in range(2):
                for j in range(2):
                    tmp[i][j] = jacobian[k][t]
                    t += 1

            xD = np.linalg.inv(tmp)
            inv_jacobian[k][0] = xD[0][0]
            inv_jacobian[k][1] = xD[0][1]
            inv_jacobian[k][2] = xD[1][0]
            inv_jacobian[k][3] = xD[1][1]
            t = 0

        return jacobian, dN_dKsi, dN_dEta, determinant, inv_jacobian

    def H_matrix(self):
        dN_dY = np.zeros((4, self.integration_points_count), float)
        dN_dX = np.zeros((4, self.integration_points_count), float)

        # Wiersz - funkcja kształtu
        # Koluma - Punkt całkowania
        # W jednym wierszu mamy wartości jednej funkcji kształtu we wszystkich punktach całkowania
        for i in range(4):
            for k in range(self.integration_points_count):
                dN_dX[i][k] = (self.dN_dKsi[i][k] * self.inv_jac[0][0] + self.dN_dEta[i][k] * (self.inv_jac[0][1]))
                dN_dY[i][k] = (self.dN_dKsi[i][k] * self.inv_jac[0][2] + self.dN_dEta[i][k] * (self.inv_jac[0][3]))

        H = []
        for i in range(self.integration_points_count):
            tmp1 = np.outer(dN_dX[:, i], np.transpose(dN_dX[:, i])) * self.determinant[i]
            tmp2 = np.outer(dN_dY[:, i], np.transpose(dN_dY[:, i])) * self.determinant[i]
            H.append(global_data.k * (tmp1 + tmp2))

        H_almost_end = self.integral(H)
        H_end = 0
        for i in range(self.integration_points_count):
            H_end += H_almost_end[i]

        return H_almost_end, H_end

    def C_matrix(self):

        N = np.zeros((4, self.integration_points_count), float)
        for i in range(int(self.integration_points_count)):
            for j in range(4):
                if j == 0:
                    N[j][i] = 0.25 * (1 - self.ksi[i]) * (1 - self.eta[i])
                elif j == 1:
                    N[j][i] = 0.25 * (1 + self.ksi[i]) * (1 - self.eta[i])
                elif j == 2:
                    N[j][i] = 0.25 * (1 + self.ksi[i]) * (1 + self.eta[i])
                elif j == 3:
                    N[j][i] = 0.25 * (1 - self.ksi[i]) * (1 + self.eta[i])

        C = []
        for i in range(self.integration_points_count):
            tmp1 = np.outer(N[:, i], np.transpose(N[:, i])) * self.determinant[i]
            C.append(global_data.c * global_data.ro * tmp1)

        C_almost_end = self.integral(C)
        C_end = 0
        for i in range(self.integration_points_count):
            C_end += C_almost_end[i]

        return C_almost_end, C_end
