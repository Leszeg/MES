from lab_1 import Plot, nN, nodes
from math import sqrt
from lab_1.Node import Node

# Drukuje wykres siatki
# Plot.plot(nN, nodes)

print('Wybierz funkcję:')
# Trzy punkty 1D
print('1. 2x^2 + 0.1x + 3')  # na zajęciach wyszło 7.3175 dla 3 punktów

# Dwa punkty 2D
print('2. -2yx^2 + 2xy + 4')  # na zajęciach wyszło 16 dla 2 punktów

# Trzy punkty 2D
print('3. -5yx^2 + 2xy^2 + 10')  # na zajęciach wyszło 40 dla 3 punktów


def integral(fun_numb):
    result = 0
    fpc = []
    nodes = []
    tmp = 0

    if int(fun_numb) == 1:
        node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
        Ak = [5 / 9, 8 / 9, 5 / 9]
        for i in range(3):
            nodes.append(Node(node_v[i], node_v[i]))
            fpc.append(2 * nodes[i].x ** 2 + 0.1 * nodes[i].x + 3)
            result = result + fpc[i] * Ak[i]
            tmp += 1
        return result

    elif int(fun_numb) == 2:

        node_v = [-1 / sqrt(3), 1 / sqrt(3)]
        Ak = [1, 1]
        for j in range(2):
            for i in range(2):
                nodes.append(Node(node_v[i], node_v[j]))
                fpc.append(-2 * nodes[tmp].y * nodes[tmp].x ** 2 + 2 * nodes[tmp].y * nodes[tmp].x + 4)
                result = result + fpc[tmp] * Ak[i] * Ak[j]
                tmp += 1
        return result


    elif int(fun_numb) == 3:
        node_v = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
        Ak = [5 / 9, 8 / 9, 5 / 9]
        for j in range(3):
            for i in range(3):
                nodes.append(Node(node_v[i], node_v[j]))
                fpc.append(-5 * nodes[tmp].y * nodes[tmp].x ** 2 + 2 * nodes[tmp].y ** 2 * nodes[tmp].x + 10)
                result = result + fpc[tmp] * Ak[i] * Ak[j]
                tmp += 1
        return result


for i in range(1, 4):
    print(integral(i))
