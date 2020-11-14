from MES import Node as n, GlobalData as g, Read as r, Plot as p, FEM_Grid as f, Element as e, Plot


# Labolatoria 1
def lab_1():
    # Czytam dane potrzebne do wygenerowania siatki z pliku
    data = r.read()

    # Laduje do klasy GlobalData
    global_data = g.GlobalData(data["W"], data["H"], data["nH"], data["nW"])

    # Odległości między węzłami na odpowiednich osiach
    d_x = global_data.W / (global_data.nW - 1)
    d_y = global_data.H / (global_data.nH - 1)

    # Inicjowanie pustych tablic
    nodes = []
    elements = []

    # Tworzenie współrzędnych węzłów (dwie pierwsze zerowe współrzędne to ksi i eta)
    for i1 in range(global_data.nW):
        for j1 in range(global_data.nH):
            nodes.append(n.Node(0, 0, round(i1 * d_x, 4), round(j1 * d_y, 4)))

    helpp = 0
    j = 0
    # Tworzenie elementów
    # Pętla for ma dodatek '+ int((global_data.nW) / 2)' ponieważ przy tworzeniu siatki
    # elementy muszą być odpowiednio numerowane i będą dodatkowe 'puste przebiegi'
    # aby zachować odpowiednią numeracja przy końcach i początkach kolumn siatki
    for i in range(global_data.nE + int((global_data.nW) / 2)):
        if helpp < global_data.nH - 1:
            ID = []
            ID.append(j)
            ID.append(ID[0] + global_data.nH)
            ID.append(ID[1] + 1)
            ID.append(ID[0] + 1)
            elements.append(e.Element(4, ID))
            helpp += 1
            j += 1
        else:
            j += 1
            helpp = 0

    # Generuje siatkę złożoną z węzłów i elementów
    fem_grid = f.FEM_Grid(nodes, elements)

    for i in range(global_data.nN):
        print(fem_grid.ND[i].x, " ", fem_grid.ND[i].y)

    for i in range(global_data.nE):
        print(f"Element {i}")
        for j in range(4):
            print(fem_grid.ELEM[i].nodes_ID[j])

    # Drukuje wykres siatki
    Plot.plot(global_data.nN, nodes)
# Opisy klas:
# FEM_Grid - zawiera listę elementów oraz węzłów
#     Element nie zna współrzędnych swoich węzłów ale zna ich ID
#     Węzły nie wiedzą do jakiego elementu należą i nie znają swojego ID
#     Klasa FEM_Grid łączy ID węzłów z elementu z ich wartościami
#     Jeżeli chcemy poznać współrzędne węzła z danego elementu musimy to zrobić poprzez klasę FEM_Grid
# GlobalData - dane odczytane z pliku
# Element - Zawiera
#     tablice z ID swoich węzłów(nodes_ID)
# Node - zawiera współrzędne punktów
# Plot - drukuje wykres FEM_Grid
# Read - czyta dane z pliku


# Labolatoria 2
def lab_2():
    element1 = e.Element(4)
    element2 = e.Element(9)
    print(element1.integral())
    print(element2.integral())


def lab_3():
    element = e.Element(4)
    element.jacobian()


def lab_4():
    element = e.Element(4)
    element.H_matrix()


def lab_5():
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

    # for i in range(nN):
    #     print(fem_grid.n[i].x, " ", fem_grid.n[i].y)
    #
    # for i in range(nE):
    #     print(f"Element {i}")
    #     for j in range(4):
    #         print(fem_grid.e[i].nodes_ID[j])

    # Drukuje wykres siatki
    Plot.plot(nN, nodes)
    for i in range(nE + int((global_data.nW) / 2)):
        elements[i].H_matrix()


lab_1()
# lab_2()
# lab_3()
lab_4()
# lab_5()
