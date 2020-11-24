from MES import Read, GlobalData as g

# Czytam dane potrzebne do wygenerowania siatki z pliku
data = Read.read()

# Laduje do klasy GlobalData
global_data = g.GlobalData(data["W"], data["H"], data["nH"], data["nW"], data["k"], data["ro"], data["c"], data["t0"])

