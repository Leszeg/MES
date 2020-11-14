# Przechowuje listę elementów i listę węzłów potrzebne do stworzenia siatki
class FEM_Grid:

    def __init__(self, N, E):
        self.ND = list(N)
        self.ELEM = list(E)

