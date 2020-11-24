# Węzeł zawiera swoje współrzędne w układach globalnym i lokalnym
class Node:
    def __init__(self, X, Y, T0):
        self.x = float(X)
        self.y = float(Y)
        self.t0 = float(T0)