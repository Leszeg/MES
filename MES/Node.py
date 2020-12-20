class Node:
    """
    Class represents the node in the global coordinate system:


    Attributes:
    ----------
    x : float
        x coordinate.
    y : float
        y coordinate
    t0 : float
        Node temperature
    bc : float
        Flag needed to check if there is a boundary condition
    """

    def __init__(self, X: float, Y: float, T0: float, BC: bool):
        """
        Constructs all the necessary attributes for the Node object

        Parameters
        ----------
        X : float
            x coordinate
        Y : float
            y coordinate
        T0 : float
            Node temperature
        BC : bool
            Flag needed to check if there is a boundary condition
        """
        self.x = X
        self.y = Y
        self.t0 = T0
        self.bc = BC
