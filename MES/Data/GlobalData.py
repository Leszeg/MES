class GlobalData:

    def __init__(self, initial_temperature, simulation_time, simulation_step_time, ambient_temperature, alfa, H, B, N_H,
                 N_B, specific_heat, conductivity, density, integration_points):
        self.it = float(initial_temperature)  # Temperatura początkowa
        self.st = float(simulation_time)  # Czas symulacji
        self.sst = float(simulation_step_time)  # Czas kroku symulacji
        self.t0 = float(ambient_temperature)  # Temperatura otoczenia
        self.alfa = float(alfa)  # Współczynnik przekazywania ciepłą przez konwekcję
        self.H = float(H)  # Wysokość siatki
        self.B = float(B)  # Szerokość siatki
        self.N_H = int(N_H)  # Ilość węzłów na wysokości
        self.N_B = int(N_B)  # Ilość węzłów na szerokości
        self.nE = int((self.N_H - 1) * (self.N_B - 1))  # Ilość elementów w siatce
        self.nN = int(self.N_H * self.N_B)  # Ilość węzłów w siatce
        self.Cw = float(specific_heat)  # Ciepło właściwe
        self.k = float(conductivity)  # Przewodność cieplna
        self.ro = float(density)  # Gęstość materiału
        self.ip = int(integration_points)  # Ilość punktów całkowania dla jednego elementu
