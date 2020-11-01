import numpy as np
import matplotlib.pyplot as plt


def plot(nN, nodes):
    # Wykres siatki
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for i in range(nN):
        plt.scatter(nodes[i].x, nodes[i].y, color='blue')
        plt.annotate(i, (nodes[i].x, nodes[i].y))

    major_ticks_x = np.arange(0, 0.1, 0.0333)
    major_ticks_y = np.arange(0, 0.21, 0.05)
    ax.set_xticks(major_ticks_x)
    ax.set_yticks(major_ticks_y)
    ax.grid(which='both')
    plt.show()
