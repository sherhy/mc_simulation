import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import numpy.polynomial.polynomial as poly

# does not call it but is necessary
from mpl_toolkits.mplot3d import Axes3D


def plot_pair_correlation(x, y, name):
    font = {
        'family': 'serif',
        'color': 'darkred',
        'weight': 'normal',
        'size': 16
    }
    plt.clf()
    plt.scatter(x, y)
    plt.xlim(0, 5)

    plt.title(name, fontdict=font)
    plt.xlabel('r', fontdict=font)
    plt.ylabel('g(r)', fontdict=font)

    plt.savefig(f"./output/{name}.png")


def clear_plt():
    plt.clf()


def scatter_van_der_waals(x, y, name, ymax=0):
    font = {
        'family': 'serif',
        'color': 'darkred',
        'weight': 'normal',
        'size': 12
    }
    plt.scatter(x, y)

    plt.title(name, fontdict=font)
    plt.xlabel('v = reduced volume', fontdict=font)
    plt.ylabel('p = reduced pressure', fontdict=font)

    plt.xlim(0.8, 1.3)
    if ymax == 0:
        plt.ylim(0, int(y[0]) + 1)
    else:
        plt.ylim(0, ymax)

    plt.savefig(f"./output/{name}.png")


def animate(shelve_db, reduced_volume: str, max_cycle):
    def plot(plist):
        ax.clear()
        x = [p.pos.x for p in plist]
        y = [p.pos.y for p in plist]
        z = [p.pos.z for p in plist]
        ax.scatter(x, y, z, alpha=0.7)
        return ax

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    filename = f"{reduced_volume}_{max_cycle}.gif"
    anim = FuncAnimation(
        fig,
        lambda i: plot(shelve_db[f"{reduced_volume}_{i}"]['plist']),
        frames=np.arange(0, max_cycle),
        interval=170)
    anim.save(f"./output/{filename}", writer='imagemagick')


def poly_approx(xs, ys, deg=4):
    coefficients = poly.polyfit(xs, ys, deg)
    x = np.arange(0.85, 1.2, 0.01)
    f_fit = poly.Polynomial(coefficients)
    plt.plot(x, f_fit(x))
