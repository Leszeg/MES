from MES import Node as n, GlobalData as g, Read as r, Plot as p, FEM_Grid as f, Element as e, Plot


# Labolatoria 1
def lab_1():
    # Czytam dane z pliku i ładuje do klasy
    data = r.read()
    global_data = g.GlobalData(data["W"], data["H"], data["nH"], data["nW"])

    # Obliczam rzeczy potrzebne do obliczenia współrzędnych
    d_x = global_data.W / (global_data.nW - 1)
    d_y = global_data.H / (global_data.nH - 1)
    i = int(global_data.nW - 1)
    j = int(global_data.nH - 1)

    # Ilości węzłów i elementów oraz inicjowanie pustych tablic
    nN = global_data.nH * global_data.nW
    nE = (global_data.nH - 1) * (global_data.nW - 1)
    nodes = []
    elements = []

    # Tworzenie współrzędnych węzłów
    for i1 in range(i + 1):
        for j1 in range(j + 1):
            nodes.append(n.Node(0, 0, round(i1 * d_x, 4), round(j1 * d_y, 4)))

    # Tworzenie elementów
    helpp = 0
    help1 = 0
    j = 0
    for i in range(nE + int((global_data.nW) / 2)):
        if helpp < global_data.nH - 1:
            ID = []
            ID.append(j)
            ID.append(ID[0] + global_data.nH)
            ID.append(ID[1] + 1)
            ID.append(ID[0] + 1)
            elements.append(e.Element(4))
            elements[help1].kons(ID)
            help1 += 1
            helpp = helpp + 1
            j = j + 1
        else:
            j = j + 1
            helpp = 0

    # Generuje siatkę złożoną z węzłów i elementów
    fem_grid = f.FEM_Grid(nodes, elements)

    for i in range(nN):
        print(fem_grid.n[i].x, " ", fem_grid.n[i].y)

    for i in range(nE):
        print(f"Element {i}")
        for j in range(4):
            print(fem_grid.e[i].nodes_ID[j])

    # Drukuje wykres siatki
    Plot.plot(nN, nodes)


# Labolatoria 2
def lab_2():
    element1 = e.Element(4)
    element2 = e.Element(9)
    print(element1.integral())
    print(element2.integral())


lab_1()
lab_2()

element = e.Element(4)
element.jacobian()