from MES import Node as n, GlobalData as g, Read as r, Plot as p, FEM_Grid as f, Element as e, Plot

def main():
    # Labolatoria
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
    a = 0
    b = 5
    c = 6
    d = 1

    # Tworzy TYLKO ELEMENTY CZTERO-WĘZŁOWE
    for i in range(global_data.nE + int((global_data.nW) / 2)):
        if helpp < global_data.nH - 1:
            ID = []
            ID.append(j)
            ID.append(ID[0] + global_data.nH)
            ID.append(ID[1] + 1)
            ID.append(ID[0] + 1)
            nod = [nodes[a], nodes[b], nodes[c], nodes[d]]
            elements.append(e.Element(4, ID, nod))
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
   # Plot.plot(global_data.nN, nodes)



# Lab 2
    #print(elements[0].integral())
    no = []
    for i in range(9):
        no.append(n.Node(0, 0, 0, 0))
    no1 = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    element = e.Element(9, no1, no)
    #print(element.integral())
# Lab 3
    #elements[0].jacobian()

# Lab 4
    #elements[0].H_matrix()

# Lab 5
#     H = []
#     J = []
#     for i in range(12):
#         J.append(elements[i].jacobian())
#         H.append(elements[i].H_matrix())
    # elements[0].jacobian()
    elements[0].H_matrix()
main()

