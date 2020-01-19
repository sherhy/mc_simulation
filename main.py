import shelve
from random import seed

from monteCarlo import run_monte_carlo
from pairCorrelation import calculate_pair_correlation
from plotter import scatter_plot
from pressure import calculate_pressure


def main():
    seed('mc')
    particle_counts = [100, 125]
    temps = [1, 2, 2.74, 5, 10]

    for n in particle_counts:
        with shelve.open("mc") as mc:
            reduced_volume = f"{216 / n:.2f}"
            cycle = run_monte_carlo(mc, n=n, kt=2.74)
            name = f"{reduced_volume}_{cycle}"

            with shelve.open("rdd") as rdd:
                calculate_pair_correlation(mc, rdd, name)
                # not necessary for calculation
                scatter_plot(rdd[f"{name}_r"], rdd[f"{name}_g"], f"g(r)-{name}")

                with shelve.open("vdw") as vdw:
                    vdw[name] = [reduced_volume, calculate_pressure(mc, rdd, name)]


if __name__ == "__main__":
    main()
