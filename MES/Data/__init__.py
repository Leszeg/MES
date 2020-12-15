from pathlib import Path
from MES.Data.GlobalData import GlobalData

# Czytam dane potrzebne do wygenerowania siatki z pliku
path = Path(__file__).parent / 'data.txt'
with open(path) as file:
    data = {}
    for line in file:
        line = line.split()
        data[line[1]] = line[0]

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
                         data["conductivity"],
                         data["density"],
                         data["integration_points"])


def print_grid_data():
    print("GLOBAL DATA")
    attr = vars(global_data)
    for keys, values in attr.items():
        print(f"{keys}" + "     " + f"{values}")
