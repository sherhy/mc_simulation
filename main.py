import shelve
from random import seed

from monteCarlo import run_monte_carlo
from pairCorrelation import calculate_pair_correlation
from plotter import plot_pair_correlation
from pressure import calculate_pressure
from van_der_waals import plot_van_der_waals


def main():
    seed('mc')
    particle_counts = [180, 186, 189, 192, 196, 199, 202, 206, 210, 213, 216, 221, 227, 230, 233, 236, 240, 243, 246, 254, ]
    temps = [0.25, 0.50, 1.00, 2.74]

    for kt in temps:
        for n in particle_counts:
            print(f"calculating for kt: {kt} n: {n}")
            with shelve.open(f"./db/mc_{kt}") as mc:
                reduced_volume = f"{216 / n:.2f}"
                with shelve.open(f"./db/cycles_{kt}") as cyc:
                    cycle = cyc[f"{reduced_volume}"] if reduced_volume in list(cyc.keys()) else 0

                cycle = run_monte_carlo(mc, n=n, kt=kt) if cycle == 0 else cycle
                name = f"{reduced_volume}_{cycle}"

                with shelve.open(f"./db/rdd_{kt}") as rdd:
                    if f"{name}_g" not in list(rdd.keys()):
                        calculate_pair_correlation(mc, rdd, name)

                    # not necessary for calculation
                    graph_name = f"Pair correlation function v={reduced_volume}"
                    plot_pair_correlation(rdd[f"{name}_r"], rdd[f"{name}_g"], graph_name)

                    with shelve.open(f"./db/vdw_{kt}") as vdw:
                        pressure = calculate_pressure(mc, rdd, name) if name not in list(vdw.keys()) else vdw[name][1]
                        print(f"pressure: {pressure:.2f}")
                        vdw[name] = [float(reduced_volume), pressure]

    plot_van_der_waals(temps)


if __name__ == "__main__":
    main()
    # tps = [0.50, 1.00, 2.74]
    # plot_van_der_waals(tps)
