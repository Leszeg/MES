from typing import List
from MES.Data import global_data
from math import sqrt
import numpy as np
from MES.Node import Node


class Element:
    """
    The class represents an element in the FEM mesh
    Creating an element involves the calculation of all its local parameters required for further global aggregations

    Attributes
    ----------
        nodes_ID : List[int]
            Node ID list in the element (in set order)

        nodes : list[Node]
            List of node coordinates in element (in set order)

        Ak : List[float]
            List of parameters needed for integration

        x : List[float]
            List of x coordinates for each node in element (in set order)

        y : List[float]
            List of y coordinates for each node in element (in set order)

        ksi : List[float]
            List of ksi coordinates for each integration_points in element (in row order)

        eta : List[float]
            List of eta coordinates for each integration_points in element (in row order)

        jacob_matrix : ndarray
            Matrix where each row is Jacobi matrix.
            To visualize inv_jac, we take elements from the row representing the Jacobi matrix
            Then insert them one by one with rows from the top

        dN_dEta : ndarray
            Matrix of dN / dEta derivatives.
            Row - Shape function
            Column - integration point
            In one line we have the values of one shape function at all integration points

        dN_dKsi : ndarray
            Matrix of dN / dKsi derivatives.
            Row - Shape function
            Column - integration point
            In one line we have the values of one shape function at all integration points

        determinant : List[float]
            List of Jacobi matrices determinant

        inv_jac : List[float]
            Matrix where each row is inverted Jacobi matrix.
            To visualize inv_jac, we take elements from the row representing the inverted matrix
            Then insert them one by one with rows from the top

        H_matrix_for_ip : list of ndarray
            List of H_matrix for each integration point in element

        H_matrix_for_element : ndarray
            H_matrix for element

        C_matrix_for_ip : ndarray
            C_matrix for each integration point in element

        C_matrix_for_element : ndarray
            C_matrix for element

        P_matrix_for_element : ndarray
            P_matrix for element

        BCH_matrix_for_element : ndarray
            BCH_matrix for element
    """

    def __init__(self, Id: List[int], nodes_c: List[Node]):
        """
        Constructs all the necessary attributes for the Element object

        Parameters
        ----------
        Id : list of int
            Node ID list in the element (in set order)

        nodes_c : list of Node
            List of node coordinates in element (in set order)
        """
        self.nodes_ID = Id

        self.nodes = nodes_c
        self.x = []
        self.y = []
        for i in self.nodes:
            self.x.append(i.x)
            self.y.append(i.y)

        if global_data.ip == 4:
            self.Ak = [1, 1]

        elif int(global_data.ip) == 9:
            self.Ak = [5 / 9, 8 / 9, 5 / 9]

        elif int(global_data.ip) == 16:
            self.Ak = [0.347855, 0.652145, 0.652145, 0.347855]

        # Współrzędne punktów całkowania w zależności od schematu
        if global_data.ip == 4:
            node_v = [-1 / sqrt(3), 1 / sqrt(3)]

        elif global_data.ip == 9:
            node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]

        elif global_data.ip == 16:
            node_v = [-0.861136, -0.339981, 0.339981, 0.861136]
        else:
            raise ValueError("Wrong value of global_data.ip")

        # Tablica współrzędnych punktów całkowania - tworzona wierszami dla elementu
        k = 0
        t = 0
        self.ksi = []
        self.eta = []
        for i in range(global_data.ip):
            self.ksi.append(node_v[k])
            self.eta.append(node_v[t])
            k += 1
            if k % len(node_v) == 0:
                t += 1
                k = 0
        self.jacob_matrix, self.dN_dEta, self.dN_dKsi, self.determinant, self.inv_jac = self.jacobian()
        self.H_matrix_for_ip, self.H_matrix_for_element = self.H_matrix()
        self.C_matrix_for_ip, self.C_matrix_for_element = self.C_matrix()
        self.P_matrix_for_element, self.BCH_matrix_for_element = self.boundary_condition()
        self.H_matrix_for_element += self.BCH_matrix_for_element

    def integral(self, matrix):
        """
        Function multiply H or C matrices (must be square) with parameters is set order

        Parameters
        ----------
        matrix : ndarray
            Square matrix H or C

        Returns
        -------
        list
            Returns the list of matrix after multiplication by parameters
        """
        result = []
        k = 0
        for j in self.Ak:
            for i in self.Ak:
                result.append(matrix[k] * i * j)
                k += 1
        return result

    def jacobian(self):
        """
        Function calculate jacobian and required parameters

        Returns
        -------
        tuple
            Return tuple of jacobian, dN_dKsi, dN_dEta, determinant and inverted jacobian
        """
        # Obliczam pochodne funkcji kształtu po ksi i eta
        dN_dKsi = np.zeros((4, global_data.ip), float)
        dN_dEta = np.zeros((4, global_data.ip), float)
        for i in range(0, global_data.ip):
            for j in range(0, global_data.ip):
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

        jacobian = np.zeros((global_data.ip, 4), float)
        jac_tmp = np.zeros(4, float)
        determinant = []
        inv_jacobian = np.zeros((global_data.ip, 4), float)

        # Cztery jakobiany ułożone wierszami
        for i in range(0, global_data.ip):
            jac_tmp[0] = np.sum(dN_dKsi[:, i] * self.x)
            jac_tmp[1] = np.sum(dN_dKsi[:, i] * self.y)
            jac_tmp[2] = np.sum(dN_dEta[:, i] * self.x)
            jac_tmp[3] = np.sum(dN_dEta[:, i] * self.y)
            jacobian = np.vstack((jacobian[:i], jac_tmp))

        for i in range(global_data.ip):
            determinant.append(jacobian[i][0] * jacobian[i][3] - jacobian[i][1] * jacobian[i][2])

        # zamiana wiersza jakobiany na macierz kwadratową
        tmp_jacobian = []
        for i in range(global_data.ip):
            tmp_jacobian.append(np.zeros((2, 2), float))

        # Obliczanie odwrotnego jakobianu
        # zamiana wiersz na macierz 2x2 i odwrócenie
        t = 0
        for k in range(global_data.ip):
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
        """
        Function calculate H_matrix

        Returns
        -------
        tuple
            Return tuple of H_matrix for each integration point and H_matrix for element
        """

        dN_dX = (self.dN_dKsi * self.inv_jac[0][0] + self.dN_dEta * (self.inv_jac[0][1]))
        dN_dY = (self.dN_dKsi * self.inv_jac[0][2] + self.dN_dEta * (self.inv_jac[0][3]))

        H = []
        for i in range(global_data.ip):
            tmp1 = np.outer(dN_dX[:, i], np.transpose(dN_dX[:, i])) * self.determinant[i]
            tmp2 = np.outer(dN_dY[:, i], np.transpose(dN_dY[:, i])) * self.determinant[i]
            H.append(global_data.k * (tmp1 + tmp2))

        H_almost_end = self.integral(H)
        H_end = 0
        for i in H_almost_end:
            H_end += i

        return H_almost_end, H_end

    @staticmethod
    def _shape_func_value(ksi, eta, Flag):
        """
            Function calculate values of shape function
            Row - shape function
            Column - iteration point
        Parameters
        ----------
        ksi : List[floats]
        eta : List[floats]
        Flag : int

        Returns
        -------
        list
            Returns matrix with values of shape function
        """

        if Flag != 0:
            r = Flag
            N = np.zeros((4, r), float)
        else:
            r = global_data.ip
            N = np.zeros((4, global_data.ip), float)

        for i in range(r):
            N[0, i] = 0.25 * (1 - ksi[i]) * (1 - eta[i])
            N[1, i] = 0.25 * (1 + ksi[i]) * (1 - eta[i])
            N[2, i] = 0.25 * (1 + ksi[i]) * (1 + eta[i])
            N[3, i] = 0.25 * (1 - ksi[i]) * (1 + eta[i])
        return N

    def C_matrix(self):
        """
        Function calculate C_matrix

        Returns
        -------
        tuple
            Return tuple of C_matrix for each integration point and C_matrix for element
        """

        N = self._shape_func_value(self.ksi, self.eta, 0)
        C = []
        for i in range(global_data.ip):
            tmp1 = np.outer(N[:, i], np.transpose(N[:, i])) * self.determinant[i]
            C.append(global_data.Cw * global_data.ro * tmp1)

        C_almost_end = self.integral(C)
        C_end = 0
        for i in C_almost_end:
            C_end += i

        return C_almost_end, C_end

    def boundary_condition(self):
        """
        Function calculate boundary_condition in element

        Returns
        -------
        tuple
            Return tuple of P_matrix and BCH_matrix
        """
        if global_data.ip == 4:
            choice2 = 2
            ksi = [-1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3), -1, - 1]
            eta = [-1, -1, -1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3)]
            N = self._shape_func_value(ksi, eta, 8)
        elif global_data.ip == 9:
            choice2 = 3
            ksi = [-sqrt(3 / 5), 0, sqrt(3 / 5), 1, 1, 1, sqrt(3 / 5), 0, -sqrt(3 / 5), -1, -1, -1]
            eta = [-1, -1, -1, -sqrt(3 / 5), 0, sqrt(3 / 5), 1, 1, 1, sqrt(3 / 5), 0, -sqrt(3 / 5)]
            N = self._shape_func_value(ksi, eta, 12)
        elif global_data.ip == 16:
            choice2 = 4
            ksi = [-0.861136, -0.339981, 0.339981, 0.861136, 1, 1, 1, 1, 0.861136, 0.339981, -0.339981, -0.861136,
                   -1, -1, -1, -1]
            eta = [-1, -1, -1, -1, -0.861136, -0.339981, 0.339981, 0.861136, 1, 1, 1, 1, 0.861136, 0.339981,
                   -0.339981, -0.861136]
            N = self._shape_func_value(ksi, eta, 16)
        else:
            raise ValueError("Wrong value of global_data.ip")

        BC = []
        Pl = []
        tmp = 0
        tmp2 = 0
        k = 0
        L_x = global_data.B / (global_data.N_B - 1)
        L_y = global_data.H / (global_data.N_H - 1)
        for i in range(3):
            if self.nodes[i].bc == 1 and self.nodes[i + 1].bc == 1:
                for j in range(choice2):
                    tmp2 += N[:, k] * self.Ak[j]
                    tmp += np.outer(N[:, k], np.transpose(N[:, k])) * self.Ak[j]

                    k += 1
                if i == 0 or i == 2:
                    BC.append(global_data.alfa * tmp * (L_x / 2))
                    Pl.append(-global_data.alfa * 1200 * tmp2 * (L_x / 2))
                else:
                    BC.append(global_data.alfa * tmp * (L_y / 2))
                    Pl.append(-global_data.alfa * 1200 * tmp2 * (L_y / 2))
                tmp = 0
                tmp2 = 0

            else:
                k += choice2

            if i == 2 and self.nodes[3].bc == 1 and self.nodes[0].bc == 1:
                for j in range(choice2):
                    tmp2 += N[:, k] * self.Ak[j]
                    tmp += np.outer(N[:, k], np.transpose(N[:, k])) * self.Ak[j]
                    k += 1
                BC.append(global_data.alfa * tmp * (L_y / 2))
                Pl.append(-global_data.alfa * 1200 * tmp2 * (L_y / 2))
                tmp = 0
                tmp2 = 0

        BCH = np.zeros((4, 4), float)
        P = np.zeros(4, float)
        for i in range(len(Pl)):
            BCH += BC[i]
            P += Pl[i]

        return P, BCH
