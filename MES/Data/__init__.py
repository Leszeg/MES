"""
The file initializes global data
"""
from pathlib import Path
from MES.Data.GlobalData import GlobalData

# Czytam dane potrzebne do wygenerowania siatki z pliku
path = Path(__file__).parent / 'data.txt'
with open(path) as file:
    data = {}
    for line in file:
        line = line.split()
        data[line[1]] = line[0]
try:
    if data["integration_points"] != '4' and data["integration_points"] != '9' and data["integration_points"] != '16':
        raise ValueError(" Ilość punktów całkowania musi wynosić 4, 9 lub 16")
except ValueError:
    print("\n Ilość punktów całkowania musi wynosić 4, 9 lub 16")
    print(" Popraw ostatni wiersz w pliku data.txt")
    exit()

# Ładuje do klasy GlobalData
global_data = GlobalData(data["initial_temperature"],
                         data["simulation_time"],
                         data["simulation_step_time"],
                         data["ambient_temperature"],
                         data["alfa"],
                         data["H"],
                         data["B"],
                         data["N_H"],
                         data["N_B"],
                         data["specific_heat"],
                         data["thermal_paste_specific_heat"],
                         data["conductivity"],
                         data["thermal_paste_conductivity"],
                         data["density"],
                         data["thermal_paste_density"],
                         data["integration_points"],
                         data["number_of_nodes_per_mm"])
