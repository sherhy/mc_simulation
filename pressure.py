import math
import shelve

from numpy import arange


def calculate_pressure(mc, rdd, name):
    def ljp(r):
        return 4 * ((1 / r) ** 12 - (1 / r) ** 6)

    rho = mc[name]['n'] / (int(mc[name]['border']) * 2) ** 3
    kt = mc[name]['dlessTemp']

    gs = rdd[name + "_g"]
    rs = rdd[name + "_r"]
    increments = [rs[i] - rs[i - 1] for i in range(1, len(rs))]

    integral = sum([rs[i] ** 3 * ljp(rs[i]) * gs[i] for i in range(0, len(rs) - 1)])

    res = rho * kt - math.pi * (2 / 3) * rho ** 2 * integral
    print(res)
    return res


if __name__ == "__main__":
    df = shelve.open('mc')
    rad_cor = shelve.open('rdd')
    vdw = shelve.open('vanDerWaals')

    df.close()
    rad_cor.close()
    vdw.close()
