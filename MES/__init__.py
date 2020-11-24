from MES import Read, GlobalData as g

# Czytam dane potrzebne do wygenerowania siatki z pliku
data = Read.read()

# Laduje do klasy GlobalData
global_data = g.GlobalData(data["W"], data["H"], data["nH"], data["nW"], data["k"], data["ro"], data["c"], data["t0"])

from MES import Node as n, FEM_Grid as f, Element as e

# Odległości między węzłami na odpowiednich osiach
d_x = global_data.W / (global_data.nW - 1)
d_y = global_data.H / (global_data.nH - 1)

# Inicjowanie pustych tablic
nodes = []
elements = []

# Tworzenie współrzędnych węzłów (dwie pierwsze zerowe współrzędne to ksi i eta)
for i1 in range(global_data.nW):
    for j1 in range(global_data.nH):
        nodes.append(n.Node(i1 * d_x, j1 * d_y, global_data.t0))

# Tworzenie elementów
# Pętla for ma dodatek '+ int((global_data.nW) / 2)' ponieważ przy tworzeniu siatki
# elementy muszą być odpowiednio numerowane i będą dodatkowe 'puste przebiegi'
# aby zachować odpowiednią numeracja przy końcach i początkach kolumn siatki
a = 0
b = int(global_data.nH)
c = int(global_data.nH + 1)
d = 1
helpp = 0
j = 0
for i in range(global_data.nE + int((global_data.nW) / 2)):
    if helpp < global_data.nH - 1:
        ID = []
        ID.append(j)
        ID.append(ID[0] + global_data.nH)
        ID.append(ID[1] + 1)
        ID.append(ID[0] + 1)
        nod = [nodes[a], nodes[b], nodes[c], nodes[d]]
        elements.append(e.Element(9, ID, nod))
        helpp += 1
        j += 1
    else:
        j += 1
        helpp = 0

grid = f.FEM_Grid(nodes, elements)
