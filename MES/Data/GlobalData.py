class GlobalData:
    """
    Class contains the global data needed for simulation

    Attributes
    ----------
        it : float
            initial_temperature

        st : float
            simulation_time

        sst : float
            simulation_step_time

        t0 : float
            ambient_temperature

        alfa : float
            convective heat transfer coefficient

        H : float
            Grid height

        B : float
            Grid Width

        N_H : int
            Number of nodes in height

        N_B : int
            Number of nodes in width

        npm : int
            number of nodes per mm

        nN : int
            number of nodes in the mesh

        nE : int
            number of elements in the mesh

        Cw : float
            Specific heat of copper

        tpCw : float
            Specific heat of thermal paste

        k : float
            Copper's thermal conductivity

        tpk : float
            Thermal paste thermal conductivity

        ro : float
            Copper density

        tpro : float
            Thermal paste density

        ip : int
            number of integration points in element
    """

    def __init__(self, initial_temperature, simulation_time, simulation_step_time,
                 ambient_temperature, alfa, H, B, N_H,
                 N_B, specific_heat, thermal_paste_specific_heat, conductivity,
                 thermal_paste_conductivity, density, thermal_paste_density,
                 integration_points, np):
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
        self.nN = int(666)  # Ilość węzłów w siatce
        self.nE = 568  # Ilość elementów w siatce
        self.Cw = float(specific_heat)  # Ciepło właściwe
        self.tpCw = float(thermal_paste_specific_heat)  # Ciepło właściwe
        self.k = float(conductivity)  # Przewodność cieplna
        self.tpk = float(thermal_paste_conductivity)  # Przewodność cieplna
        self.ro = float(density)  # Gęstość materiału
        self.tpro = float(thermal_paste_density)  # Gęstość materiału
        self.ip = int(integration_points)  # Ilość punktów całkowania dla jednego elementu
