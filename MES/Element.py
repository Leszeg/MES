from math import sqrt
from MES.Node import Node


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
