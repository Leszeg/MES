# Pakiet MES

## Rozwiązuje równanie Fouriera nieustalonego przepływu ciepła dla elementu prostokątnego za pomocą metody elementów skończonych

## Założenia:

- Badany materiał jest izotropowy
- badamy siatkę 2D

## Wywołanie

Do pliku main.py importujemy pakiet MES i wywołujemy metodę run()
Ze względu na użycie wielowątkowości plik musi zawierać klauzulę sprawdzającą czy metoda wywołana jest bezpośrednio z
pliku main Przykładowy plik main.py

```python
import MES

if __name__ == '__main__':
    MES.run()
```

## Dane

Wewnątrz pakietu MES znajduję się pakie Data a w nim plik data.txt. Tutaj można wprowadzić własne parametry symulacji
podmieniając wartości liczbowe

```text
100     initial_temperature
500     simulation_time
50      simulation_step_time
1200    ambient_temperature
300     alfa
0.100   H
0.100   B
4       N_H
4       N_B
700     specific_heat
25      conductivity
7800    density
4       integration_points
```

## Działanie programu:

- Tworze globalną siatkę składającą się z węzłów
- Tworze w siatce elementy przypisując im odpowiednie węzły
    - Węzły numeruje w odpowiedni sposób
    - Każdy element zależnie od wybranego schematu(4,9,16) ustala współrzędne lokalne punktów całkowania w układzie
      ksi-eta
    - Obliczam odwrócony jakobian aby uzależnić funkcje kształtu od zmiennych globalnych a nie lokalnych
    - Obliczam macierze H i C dla każdego elementu
        - H - całka po objętości przemnożona przez współczynnik przewodzenia ciepła
            - składnik 1 - wektor pochodnych funkcji kształtu po x przemnożony przez ten sam wektor ale transponowany
            - składnik 2 - wektor pochodnych funkcji kształtu po y przemnożony przez ten sam wektor ale transponowany
            - macierz H dla punktu całkowania - dodanie składnika 1 i składnika 2 i przemnożenie przez współczynnik
              przewodzenia ciepłą
            ```python
                    for i in range(global_data.ip):
                        tmp1 = np.outer(dN_dX[:, i], np.transpose(dN_dX[:, i])) * self.determinant[i]
                        tmp2 = np.outer(dN_dY[:, i], np.transpose(dN_dY[:, i])) * self.determinant[i]
                        H.append(global_data.k * (tmp1 + tmp2))
            
                    H_almost_end = self.integral(H)
                    H_end = 0
                    for i in H_almost_end:
                        H_end += i
            
                    return H_almost_end, H_end
            ```
        - C - całka po objętości przemnożona przez gęstość i pojemność cieplną
            - pod całką znajduje się wektor wartości funkcji kształtu przemnożony przez ten sam wektor ale transponowany
            ```python
                    for i in range(global_data.ip):
                        tmp1 = np.outer(N[:, i], np.transpose(N[:, i])) * self.determinant[i]
                        C.append(global_data.Cw * global_data.ro * tmp1)
            
                    C_almost_end = self.integral(C)
                    C_end = 0
                    for i in C_almost_end:
                        C_end += i
                    return C_almost_end, C_end
            ```
    - całkowanie zrealizowane poprzez pomnożenie przez wyznacznik macierzy Jacobiego i przemnożenie przez odpowiednie
      wagi (kwadratury Gaussa)
    - Warunek brzegowy (Dirichleta - konwekcyjny) rozbijam na macierz Hbc i wektor P
        - Obliczam macierz Hbc i dodaje ją do macierzy H
        - Obliczam wektor obciążeń P który wykorzystam później
- Po obliczeniu wszystkich elementów agreguje ich macierze lokalne (H, C, P) do macierzy globalnych dla całej siatki
- Teraz obliczenie temperatury w węzłach sprowadza sie do rozwiązania układu macierzowego
    - Obliczam Hz, Pz
        - Hz - globalna macierz H do której dodajemy globalną macierz C podzieloną przez krok czasowy symulacji
        - Pz - składnik 3 dodajemy do macierzy P i całość mnożymy razy -1
            - mnożę przez -1 bo przenoszę P na prawą stronę równania
            - składnik 3 - globalna macierz C podzielona przez krok czasowy symulacji a następnie wszystko przemnożone
              przez wektor temperatur
- Rozwiązuje układ Hz * {t} = Pz
- Obliczony wektor temperatur podstawiam do węzłów w siatce
- Symulacja polega na ponownym policzeniu macierzy H;C;P, ich agregacji, rozwiązaniu układu i podstawieniu temperatur do
  węzłów

```python
    for i in range(0, int(global_data.st) + int(global_data.sst), int(global_data.sst)):
    global g
if i == 0:
    g.get_temps(True, temps)
    temps = g.temperature_of_nodes
    g.to_file()
    pool.map(recalculate_matrix, range(len(g.elements)))
else:
    pool.map(recalculate_matrix, range(len(g.elements)))
g.matrix_aggregation()
d[f'Time {i}'] = temps[0]

Hz = g.H_global + (g.C_global / global_data.sst)
x = -dot(g.C_global / global_data.sst, temps.T).T
Pz = -(g.P_global + x)
temps = solve(Hz, Pz.T).T
```

# Pliki wynikowe

Program zwraca dwa pliki

- MES.xls - Arkusz excel zawieracjący informacje o każdym elemencie przed rozpoczęciem symulacji (wartości początkowe
  siatki)
    - Plik ma na celu sprawdzenie poprawności obliczeń
- Temperatures.csv - Plik csv zawierający temperatury W każdym węźle w każdym kroku symulacji
    - Każda kolumna zawiera chwilę czasową i temperatury we wszystkich węzłach zgodnie z numeracją
    - Po otwarciu w programie Excel wszystkie dane będą w kolumnie A
    - Należy wybrać z opcję 'tekst jako kolumny' z zakładki dane i jako separator wybrać przecinek
    
