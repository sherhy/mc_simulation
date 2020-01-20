import shelve
from math import pi

from numpy import arange

from plotter import plot_pair_correlation


def get_all_distances(df):
    p_list = df['plist']
    g_list = df['ghost']
    distances = list()
    for p in p_list:
        distances += list(map(lambda point: p.pos.get_distance_to(point.pos), p_list + g_list))
    return distances


def g(n, distances, r):
    dr = 0.1
    in_range_count = len(list(filter(lambda length: abs(length - r) < dr, distances)))
    return 110 * in_range_count / (n ** 2 * 4 * pi * r ** 2 * dr)


def calculate_pair_correlation(mc, rdd, name):
    g_values = list()
    r_values = list()
    distances = get_all_distances(mc[name])
    rdd[name + "_all_dist"] = distances.copy()
    n = mc[name]['n']

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

    rdd[name + "_g"] = g_values
    rdd[name + "_r"] = r_values


if __name__ == '__main__':
    with shelve.open("mc") as db:
        with shelve.open("rdd") as radial_dist:
            reduced_volume = "1.73"
            cycle_count = "10"
            db_name = reduced_volume + "_" + cycle_count

            calculate_pair_correlation(db, radial_dist, db_name)
            plot_pair_correlation(radial_dist[f'{db_name}_r'], radial_dist[f'{db_name}_g'], f"g(r)-{db_name}")

