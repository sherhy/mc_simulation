import shelve
from random import seed

from monteCarlo import run_monte_carlo
from pairCorrelation import calculate_pair_correlation


def main(start):
    kts = [2.74, 5, 20, 100, 1]

    end = run_monte_carlo(n=2, border_margin=1)
    calculate_pair_correlation(db, rdd, end)


if __name__ == "__main__":
    seed('mc')
    db = shelve.open("mcSimulation")
    rdd = shelve.open("radialDistribution")

    main()

    db.close()
    rdd.close()
