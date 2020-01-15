import shelve
from math import pi

from numpy import arange

from plotter import scatter_plot


def g(df, r):
    dr = 0.1
    p_list = df['plist']
    g_list = df['ghost']
    n = len(p_list)
    count = 0
    for p in p_list:
        distances = list(map(lambda point: point.pos.get_distance_to(p.pos).get_magnitude(), p_list + g_list))
        in_range = list(filter(lambda length: abs(length - r) < dr, distances))

        count += len(in_range)
    return count / (n * 4 * pi * r ** 2 * dr)


if __name__ == '__main__':
    db = shelve.open("mcSimulation")
    rdd = shelve.open("radialDistribution")

    g_values = list()
    x_values = arange(0.25, 5, 0.05)

    # change here
    cycle = '70'
    for x in x_values:
        g_val = g(db[cycle], x)
        if x.is_integer():
            print(f"r={x} ")
        g_values.append(g_val)
    rdd[cycle] = g_values

    scatter_plot(x_values, g_values, f"g(r)-{cycle}")

    db.close()
    rdd.close()
