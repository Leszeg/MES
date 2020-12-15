from MES.Data import print_grid_data, global_data
from MES.FEM_Grid import FEM_Grid
from numpy import fabs
import numpy as np


def Gauss(macierz, wektor_wartosci, n, wektor_wynikow):
    licznik = 1

    # Przekształcam macierz do postaci trójkątnej
    for k in range(0, n - 1):

        # Obsługa zer na przekątnej macierzy (zamiana wierszy)
        if fabs(macierz[k, k]) < 1.0e-12:
            for i in range(k + 1, n):
                if fabs(macierz[i, k]) > fabs(macierz[k, k]):
                    macierz[[k, i]] = macierz[[i, k]]
                    wektor_wartosci[[k, i]] = wektor_wartosci[[i, k]]
                    break
        # Koniec obsługi

        for i in range(k + 1, n):
            if macierz[i, k] == 0: continue  # Warunek pomijający 0 pod przekątną macierzy
            factor = macierz[i, k] / macierz[k, k]

            for j in range(k, n):
                macierz[i, j] = macierz[i, j] - (macierz[k, j] * factor)
            wektor_wartosci[i] = wektor_wartosci[i] - (wektor_wartosci[k] * factor)

            licznik += 1

    licznik = 1

    # Rozwiązuje układ równań "od końca"
    wektor_wynikow[n - 1] = wektor_wartosci[n - 1] / macierz[n - 1, n - 1]
    for i in range(n - 2, -1, -1):
        sum_ax = 0
        for j in range(i + 1, n):
            sum_ax += macierz[i, j] * wektor_wynikow[j]
        wektor_wynikow[i] = (wektor_wartosci[i] - sum_ax) / macierz[i, i]
        licznik += 1

    return wektor_wynikow


def get_temps(g: FEM_Grid):
    x = g.temperature_of_nodes.T
    return x.T


temps = np.zeros((global_data.N_B * global_data.N_H, 1))

for i in range(0, int(global_data.st), int(global_data.sst)):
    if i == 0:
        g = FEM_Grid(True, temps)
        temps = get_temps(g)
        # g.plot_grid()
    else:
        g = FEM_Grid(False, temps)
        temps = get_temps(g)
    Hz = g.H_global + (g.C_global / global_data.sst)
    x = -np.dot(g.C_global / global_data.sst, temps.T).T
    Pz = -(g.P_global + x)
    temps = np.linalg.solve(Hz, Pz.T).T
    p = Gauss(Hz, Pz.T, 960, np.zeros((global_data.N_B * global_data.N_H, 1)))

    print(i)
    print(max(temps[0]))
    print(min(temps[0]))
    print("\n")
