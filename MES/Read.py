# Odczytuje dane z pliku
# W   0.1 - szerokość siatki
# H   0.2 - wysokość siatki
# nH  5 - liczba węzłów na wysokości
# nW  4 - liczba węzłów na szerokości

def read():
    with open('data.txt') as file:
        tmp = {}
        for line in file:
            line = line.split()
            tmp[line[0]] = line[1]
    return tmp
