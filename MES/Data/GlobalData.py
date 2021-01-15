class GlobalData:

    def __init__(self, initial_temperature, simulation_time, simulation_step_time, ambient_temperature, alfa, H, B, N_H,
                 N_B, specific_heat, conductivity, density, integration_points, np):
        self.it = float(initial_temperature)  # Temperatura początkowa
        self.st = float(simulation_time)  # Czas symulacji
        self.sst = float(simulation_step_time)  # Czas kroku symulacji
        self.t0 = float(ambient_temperature)  # Temperatura otoczenia
        self.alfa = float(alfa)  # Współczynnik przekazywania ciepłą przez konwekcję
        self.H = float(H)  # Wysokość siatki
        self.B = float(B)  # Szerokość siatki
        self.N_H = int(N_H)  # Ilość węzłów na wysokości
        self.N_B = int(N_B)  # Ilość węzłów na szerokości
        self.npm = int(np)  # Ilość węzłów na mm
        self.nN = int(4 * (self.N_H * 2 * self.npm) + 3 * (2 * self.npm * 3))  # Ilość węzłów w siatce
        self.Cw = float(specific_heat)  # Ciepło właściwe
        self.k = float(conductivity)  # Przewodność cieplna
        self.ro = float(density)  # Gęstość materiału
        self.ip = int(integration_points)  # Ilość punktów całkowania dla jednego elementu
