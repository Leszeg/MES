class Node:
    def __init__(self, X, Y, T0, BC):
        # Współrzędne globalne - (x,y)
        self.x = float(X)
        self.y = float(Y)

        # Temperatura w węźle
        self.t0 = float(T0)

        # Określenie występowanie warunku brzegowego
        self.bc = BC
