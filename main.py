import shelve
from random import seed

from monteCarlo import run_monte_carlo
from pairCorrelation import calculate_pair_correlation
from plotter import plot_pair_correlation, plot_van_der_waals
from pressure import calculate_pressure


def main():
    seed('mc')
    particle_counts = [45, 61, 86, 120, 155, 180, 216, 221, 227, 230, 233, 236, 240, 243, 246, 254, ]
    temps = [0.50, 1.00, 2.74]

    for kt in temps:
        for n in particle_counts:
            print(f"calculating for kt: {kt} n: {n}")
            with shelve.open(f"./db/mc_{kt}") as mc:
                reduced_volume = f"{216 / n:.2f}"
                with shelve.open(f"./db/cycles_{kt}") as cyc:
                    cycle = run_monte_carlo(mc, n=n, kt=kt) \
                        if reduced_volume not in list(cyc.keys()) \
                        else cyc[f"{reduced_volume}"]
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
                        vdw[name] = [reduced_volume, pressure]

    show_van_der_waals()


def show_van_der_waals():
    print("creating Van Der Waals.png")
    temps = [1.00, 2.74]
    all_v_val = list()
    all_p_val = list()
    for kt in temps:
        v_values = list()
        p_values = list()
        with shelve.open(f"./db/vdw_{kt}") as vdw:
            skip_list = ["0.75_7"]
            keys = list(vdw.keys())
            keys.sort()
            # print(keys)
            for key in keys:
                if key in skip_list:
                    continue
                v_values.append(vdw[key][0])
                all_v_val.append(vdw[key][0])
                p_values.append(vdw[key][1])
                all_p_val.append(vdw[key][1])
        plot_van_der_waals(v_values, p_values, f"Van Der Waals Equation t={kt}")
    plot_van_der_waals(all_v_val, all_p_val, "Van Der Waals Equation")


if __name__ == "__main__":
    # main()
    show_van_der_waals()
