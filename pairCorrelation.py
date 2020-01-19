import shelve
from math import pi

from numpy import arange

from plotter import scatter_plot, poly_approx


def g(n, distances, r):
    dr = 0.1
    in_range = list(filter(lambda length: abs(length - r) < dr, distances))
    count = len(in_range)
    return count / (n * 4 * pi * r ** 2 * dr)


def get_all_distances(df):
    p_list = df['plist']
    g_list = df['ghost']
    distances = list()
    for p in p_list:
        distances += list(map(lambda point: point.pos.get_distance_to(p.pos), p_list + g_list))
    return distances


def calculate_pair_correlation(df, rdd, cycle: str):
    g_values = list()
    r_values = list()
    distances = get_all_distances(df[cycle])
    rdd[cycle + "_all_dist"] = distances
    n = df[cycle]['n']

    for r in arange(0.4, 0.8, 0.05):
        g_val = g(n, distances, r)
        g_values.append(g_val)
        r_values.append(r)

    for r in arange(0.8, 1.3, 0.01):
        g_val = g(n, distances, r)
        g_values.append(g_val)
        r_values.append(r)

    for r in arange(1.3, 5, 0.05):
        g_val = g(n, distances, r)
        g_values.append(g_val)
        r_values.append(r)

    rdd[cycle + "_g"] = g_values
    rdd[cycle + "_r"] = r_values


if __name__ == '__main__':
    db = shelve.open("mcSimulation")
    radial_dist = shelve.open("radialDistribution")

    calculate_pair_correlation(db, radial_dist, '18')

    cycle_count = 18
    # not necessary for calculation
    scatter_plot(radial_dist[f'{cycle_count}_r'], radial_dist[f'{cycle_count}_g'], f"g(r)-{cycle_count}")
    # poly_approx(r_values, g_values, f"g(r)-{cycle}-poly")

    db.close()
    radial_dist.close()
