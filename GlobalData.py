class GlobalData:
    W = 0
    H = 0
    nH = 0
    nW = 0

    def __init__(self, w, h, nh, nw):
        self.W = float(w)
        self.H = float(h)
        self.nH = int(nh)
        self.nW = int(nw)
