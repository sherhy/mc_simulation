import math
import shelve

from numpy import arange


def calculate_pressure(db, rdd, cycle: str):
    def ljp(r):
        return 4 * ((1 / r) ** 12 - (1 / r) ** 6)

    rho = 125 / 216
    kt = 2.74
    dr = .05

    rs = rdd[cycle + "_r"]
    increments = [rs[i] - rs[i-1] for i in range(1, len(rs))]
    gs = rdd[cycle + "_g"]

    integral = sum([rs[i]**3 * ljp(rs[i]) * gs[i] * (rs[i+1] - rs[i]) for i in range(0, len(rs)-1)])

    res = rho * kt - math.pi * (2/3) * rho**2 * integral
    print(res)
    return res


if __name__ == "__main__":
    df = shelve.open('mcSimulation')
    rad_cor = shelve.open('radialDistribution')
    vdw = shelve.open('vanDerWaals')

    calculate_pressure(df, rad_cor, '18')

    df.close()
    rad_cor.close()
