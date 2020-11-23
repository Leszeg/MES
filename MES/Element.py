from math import sqrt
import numpy as np
from MES.__init__ import global_data
from MES import IntegrationPoint as IP


class Element:

    def __init__(self, c, Id, nodes_c):
        self.nodes_ID = list(Id)

        # Dodane na potrzeby całkowania informacja o ilości węzłów w elemencie
        self.integration_points_count = c

        self.nodes = nodes_c

        if c == 4:
            node_v = [-1 / sqrt(3), 1 / sqrt(3)]
            self.integration_points = []

            for i in range(self.integration_points_count):
                self.integration_points.append(IP.IntegrationPoint())

            self.integration_points[0].ksi = node_v[0]
            self.integration_points[0].eta = node_v[0]

            self.integration_points[1].ksi = node_v[1]
            self.integration_points[1].eta = node_v[0]

            self.integration_points[2].ksi = node_v[1]
            self.integration_points[2].eta = node_v[1]

            self.integration_points[3].ksi = node_v[0]
            self.integration_points[3].eta = node_v[1]

        elif c == 9:
            node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
            self.integration_points = []

            for i in range(self.integration_points_count):
                self.integration_points.append(IP.IntegrationPoint())

            self.integration_points[0].ksi = node_v[0]
            self.integration_points[0].eta = node_v[0]

            self.integration_points[1].ksi = node_v[1]
            self.integration_points[1].eta = node_v[0]

            self.integration_points[2].ksi = node_v[2]
            self.integration_points[2].eta = node_v[0]

            self.integration_points[3].ksi = node_v[0]
            self.integration_points[3].eta = node_v[1]

            self.integration_points[4].ksi = node_v[1]
            self.integration_points[4].eta = node_v[1]

            self.integration_points[5].ksi = node_v[2]
            self.integration_points[5].eta = node_v[1]

            self.integration_points[6].ksi = node_v[0]
            self.integration_points[6].eta = node_v[2]

            self.integration_points[7].ksi = node_v[1]
            self.integration_points[7].eta = node_v[2]

            self.integration_points[8].ksi = node_v[2]
            self.integration_points[8].eta = node_v[2]

        elif c == 16:
            pass

    def integral(self, H_matrix):
        result = []
        fpc = []
        tmp = 0

        if self.integration_points_count == 4:
            # Wagi
            Ak = [1, 1]
            k = 0

            for i in range(2):
                for j in range(2):
                    result.append(H_matrix[k] * Ak[i] * Ak[j])
                    k += 1

            return result

        elif int(self.integration_points_count) == 9:
            # Wagi
            Ak = [5 / 9, 8 / 9, 5 / 9]
            k = 0

            for j in range(3):
                for i in range(3):
                    result.append(H_matrix[k] * Ak[i] * Ak[j])
                    k += 1
            return result

    def jacobian(self):

        # Współrzędne globalne elementu

        x = []
        y = []
        for i in range(4):
            x.append(self.nodes[i].x)
            y.append(self.nodes[i].y)

        # Współrzędne lokalne elementu numeracja jak z siatką FEM
        ksi = []
        eta = []
        for i in range(self.integration_points_count):
            ksi.append(self.integration_points[i].ksi)
            eta.append(self.integration_points[i].eta)

        # Obliczam pochodne funkcji kształtu po ksi i eta
        dN_dKsi = np.zeros((4, self.integration_points_count), float)
        dN_dEta = np.zeros((4, self.integration_points_count), float)
        for i in range(0, self.integration_points_count):
            for j in range(0, self.integration_points_count):
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

        jacobian = np.zeros((self.integration_points_count, 4), float)
        determinant = []
        inv_jacobian = np.zeros((self.integration_points_count, 4), float)

        # Cztery jakobiany 2x2 ułożone wierszami
        k = 0
        sum = 0
        for i in range(0, self.integration_points_count):
            for j in range(0, 4):
                sum += dN_dKsi[j][i] * x[j]
            jacobian[i][k] = sum
            k += 1
            sum = 0

            for j in range(0, 4):
                sum += dN_dKsi[j][i] * y[j]
            jacobian[i][k] = sum
            k += 1
            sum = 0

            for j in range(0, 4):
                sum += dN_dEta[j][i] * x[j]
            jacobian[i][k] = sum
            k += 1
            sum = 0

            for j in range(0, 4):
                sum += dN_dEta[j][i] * y[j]
            jacobian[i][k] = sum
            k = 0
            sum = 0

        for i in range(self.integration_points_count):
            determinant.append(jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2])

        tmp_jacobian = []
        for i in range(self.integration_points_count):
            tmp_jacobian.append(np.zeros((2, 2), float))
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

        # print("///////////X////////////")
        # print(x)
        # print("///////////Y////////////")
        # print(y)
        #
        # print("//////////ETA///////////")
        # for i in range(self.integration_points_count):
        #     print(eta[i])
        #
        # print("//////////KSI///////////")
        # for i in range(self.integration_points_count):
        #     print(ksi[i])
        #
        # # Pochodne w macierzy ułożone wierszami - pierwszy wiersz - dN1_dEta itd.
        # print('////////// POCHODNA FUNKCJI KSZTAŁTU ETA //////////')
        # for i in range(self.integration_points_count):
        #     print(dN_dKsi[i])
        # # Pochodne w macierzy ułożone wierszami - pierwszy wiersz - dN1_dKsi itd.
        # print('////////// POCHODNA FUNKCJI KSZTAŁTU KSI //////////')
        # for i in range(self.integration_points_count):
        #     print(dN_dEta[i])
        #
        # print('////////// JACOBIAN //////////')
        # for i in range(self.integration_points_count):
        #     print(jacobian[i])
        #
        # print('///////// DETERMINANT //////////')
        # for i in range(self.integration_points_count):
        #     print(determinant[i])
        #
        # print('////////// INVERSE JACOBIAN //////////')
        # for i in range(self.integration_points_count):
        #     print(inv_jacobian[i])

        return jacobian, dN_dKsi, dN_dEta, determinant, inv_jacobian

    def H_matrix(self):
        data = self.jacobian()
        dN_dEta = data[1]
        dN_dKsi = data[2]
        determinant = data[3]
        inv_jac = data[4]
        dN_dY = np.zeros((4, self.integration_points_count), float)
        dN_dX = np.zeros((4, self.integration_points_count), float)

        # Wiersz - funkcja kształtu
        # Koluma - Punkt całkowania
        # W jednym wierszu mamy wartości jednej funkcji kształtu we wszystkich punktach całkowania
        k = 0

        for i in range(4):
            for k in range(self.integration_points_count):
                dN_dX[i][k] = (dN_dKsi[i][k] * inv_jac[0][0] + dN_dEta[i][k] * (inv_jac[0][1]))
                dN_dY[i][k] = (dN_dKsi[i][k] * inv_jac[0][2] + dN_dEta[i][k] * (inv_jac[0][3]))

        H = []
        for i in range(self.integration_points_count):
            tmp1 = np.outer(dN_dX[:, i], np.transpose(dN_dX[:, i])) * determinant[i]
            tmp2 = np.outer(dN_dY[:, i], np.transpose(dN_dY[:, i])) * determinant[i]
            H.append(global_data.k * (tmp1 + tmp2))

        H_almost_end = self.integral(H)
        H_end = 0
        for i in range(self.integration_points_count):
            H_end += H_almost_end[i]

        print("\n\n")
        print("////////// H //////////")
        print(H_end)
        return H_end
