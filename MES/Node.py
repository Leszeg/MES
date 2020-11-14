# Węzeł zawiera swoje współrzędne w układach globalnym i lokalnym
class Node:
    def __init__(self, KSI, ETA, X, Y):
        self.ksi = float(KSI)
        self.eta = float(ETA)
        self.x = float(X)
        self.y = float(Y)
