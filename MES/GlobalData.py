class GlobalData:

    def __init__(self, w, h, nh, nw):
        self.W = float(w)
        self.H = float(h)
        self.nH = int(nh)
        self.nW = int(nw)
        self.nE = (self.nH - 1) * (self.nW - 1)
        self.nN = self.nH * self.nW
