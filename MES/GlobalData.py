class GlobalData:

    def __init__(self, w, h, nh, nw, K, RO, C, T0):
        self.W = float(w)
        self.H = float(h)
        self.nH = int(nh)
        self.nW = int(nw)
        self.nE = (self.nH - 1) * (self.nW - 1)
        self.nN = self.nH * self.nW
        self.k = float(K)
        self.ro = float(RO)
        self.c = float(C)
        self.t0 = int(T0)
