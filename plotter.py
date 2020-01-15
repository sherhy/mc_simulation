import shelve

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# does not call it but is necessary
from mpl_toolkits.mplot3d import Axes3D
from particle import Particle
from vector import Vector


def scatter_plot(x, y, name):
    plt.scatter(x, y)
    plt.savefig(f"./output/{name}.png")


def animate(db, max_cycle):
    def plot(plist):
        ax.clear()
        x = [p.pos.x for p in plist]
        y = [p.pos.y for p in plist]
        z = [p.pos.z for p in plist]
        ax.scatter(x, y, z, alpha=0.7)
        return ax

    def update(i):
        return plot(db[f"{i}"]['plist'])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    filename = f"{max_cycle}cycles.gif"
    anim = FuncAnimation(fig, update, frames=np.arange(0, max_cycle), interval=170)
    anim.save(f"./output/{filename}", writer='imagemagick')


if __name__ == '__main__':
    db = shelve.open("mcSimulation")

    # change here
    animate()

    db.close()
