import shelve
from random import seed

from monteCarlo import run_monte_carlo
from pairCorrelation import calculate_pair_correlation
from plotter import plot_pair_correlation, plot_van_der_waals
from pressure import calculate_pressure


def main():
    seed('mc')
    particle_counts = [45, 61, 86, 120, 155, 180, 216]
    temps = [1, 2, 2.74, 5, 10]

    for n in particle_counts:
        with shelve.open("mc") as mc:
            reduced_volume = f"{216 / n:.2f}"
            cycle = run_monte_carlo(mc, n=n, kt=2.74)
            name = f"{reduced_volume}_{cycle}"

            with shelve.open("rdd") as rdd:
                calculate_pair_correlation(mc, rdd, name)
                # not necessary for calculation
                graph_name = f"Pair correlation function v={reduced_volume}"
                plot_pair_correlation(rdd[f"{name}_r"], rdd[f"{name}_g"], graph_name)

                with shelve.open("vdw") as vdw:
                    vdw[name] = [reduced_volume, calculate_pressure(mc, rdd, name)]
        break


def show_van_der_waals():
    with shelve.open('vdw') as vdw:
        v_values = list()
        p_values = list()
        skip_list = ['0.72_5', '0.72_10', '1.20_7', '1.00_7', ]
        # print(list(vdw.keys()))
        keys = list(vdw.keys())
        keys.sort()
        print(keys)
        for key in keys:
            if key in skip_list:
                continue
            v_values.append(vdw[key][0])
            p_values.append(vdw[key][1])
        plot_van_der_waals(v_values, p_values, "Van Der Waals")


if __name__ == "__main__":
    main()
