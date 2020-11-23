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
        x = []
        y = []
        for i in range(4):
            x.append(self.nodes[i].x)
            y.append(self.nodes[i].y)

        # Współrzędne lokalne elementu numeracja jak z siatką FEM
        ksi = [self.integration_points[0].ksi, self.integration_points[1].ksi, self.integration_points[2].ksi, self.integration_points[3].ksi]
        eta = [self.integration_points[0].eta, self.integration_points[1].eta, self.integration_points[2].eta, self.integration_points[3].eta]

        # Obliczam pochodne funkcji kształtu po ksi i eta
        dN_dKsi = np.zeros((4, 4), float)
        dN_dEta = np.zeros((4, 4), float)
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

        jacobian = np.zeros((4, 4), float)
        determinant = []
        inv_jacobian = np.zeros((4, 4), float)
        m = dN_dEta
        dN_dEta = dN_dKsi
        dN_dKsi = m

        # Cztery jakobiany 2x2 ułożone wierszami
        for i in range(0, self.integration_points_count):
            for j in range(0, self.integration_points_count):
                for k in range(0, self.integration_points_count):
                    if j == 0:
                        jacobian[i][j] += dN_dEta[k][j] * x[k]
                    elif j == 1:
                        jacobian[i][j] += dN_dKsi[k][j] * x[k]
                    elif j == 2:
                        jacobian[i][j] += dN_dEta[k][j] * y[k]
                    elif j == 3:
                        jacobian[i][j] += dN_dKsi[k][j] * y[k]

        for i in range(self.integration_points_count):
            determinant.append(jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2])

        # tmp_jacobian = np.zeros((2, 2), float)
        #
        # tmp_jacobian[0][0] = jacobian[0][0]
        # tmp_jacobian[0][1] = jacobian[0][1]
        # tmp_jacobian[1][0] = jacobian[0][2]
        # tmp_jacobian[1][1] = jacobian[0][3]
        #
        # xD = np.linalg.inv(tmp_jacobian)
        # tmp2 = np.zeros((1, 4), float)
        #
        # tmp2[0][0] = xD[0][0]
        # tmp2[0][1] = xD[0][1]
        # tmp2[0][2] = xD[1][0]
        # tmp2[0][3] = xD[1][1]


        k = 0
        for i in range(self.integration_points_count):
            for j in range(3, -1, -1):
                if k == 0 or k == 3:
                    inv_jacobian[i][k] = round(jacobian[i][j] / determinant[i], 3)
                elif k == 1 or k == 2:
                    inv_jacobian[i][k] = round(jacobian[i][j] / determinant[i], 3)
                k += 1
            k = 0
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
        dN_dY = np.zeros((4, 4), float)
        dN_dX = np.zeros((4, 4), float)

        # Wiersz - funkcja kształtu
        # Koluma - Punkt całkowania
        # W jednym wierszu mamy wartości jednej funkcji kształtu we wszystkich punktach całkowania
        for i in range(4):
            for k in range(4):
                dN_dX[i][k] = (dN_dKsi[i][k] * inv_jac[0][0] + dN_dEta[i][k] * (inv_jac[0][1]))
                dN_dY[i][k] = (dN_dKsi[i][k] * inv_jac[0][2] + dN_dEta[i][k] * (inv_jac[0][3]))

        H = []
        for i in range(0, 4):
            tmp1 = np.outer(dN_dX[:, i], np.transpose(dN_dX[:, i])) * determinant[i]
            tmp2 = np.outer(dN_dY[:, i], np.transpose(dN_dY[:, i])) * determinant[i]
            H.append(global_data.k * (tmp1 + tmp2))

        H_almost_end = self.integral(H)

        H_end = (H_almost_end[0] + H_almost_end[1] + H_almost_end[2] + H_almost_end[3])

        print("\n\n")
        print("////////// H //////////")
        print(H_end)
        return H_end
