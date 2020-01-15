import shelve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
#does not call it but is necessary
from mpl_toolkits.mplot3d import Axes3D

def plot(plist):
    ax.clear()
    x = [p.pos.x for p in plist]
    y = [p.pos.y for p in plist]
    z = [p.pos.z for p in plist]
    ax.scatter(x, y, z, alpha=0.7)
    return ax

def update(i):
    return plot(db[f"{i}"]['plist'])

if __name__ == '__main__':
    db = shelve.open("mcSimulation")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    filename = "30cycles.gif"
    anim = FuncAnimation(fig, update, frames=np.arange(0,41), interval=170)
    # anim.save(f"./output/{filename}", writer='imagemagick')
    plt.show()

    db.close()