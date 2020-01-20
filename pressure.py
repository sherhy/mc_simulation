import math
import shelve


def calculate_pressure(mc, rdd, name):
    def ljp(r):
        return 4 * ((1 / r) ** 12 - (1 / r) ** 6)

    rho = mc[name]['n'] / (int(mc[name]['border']) * 2) ** 3
    kt = mc[name]['dlessTemp']
    sigma = 1
    r_cutoff = 5

    gs = rdd[name + "_g"]
    rs = rdd[name + "_r"]
    increments = [rs[i] - rs[i - 1] for i in range(1, len(rs))]

    integral = sum([rs[i] ** 3 * ljp(rs[i]) * gs[i] for i in range(0, len(rs) - 1)])

    pressure = rho * kt - math.pi * (2 / 3) * rho ** 2 * integral
    constant_coefficient = (16 / 9) * math.pi * rho ** 2 * sigma ** 3
    analytic_pressure = constant_coefficient * (2 * (sigma / r_cutoff) ** 9 - 3 * (sigma / r_cutoff) ** 3)
    pressure += analytic_pressure

    print(pressure)
    return pressure


if __name__ == "__main__":
    df = shelve.open('mc')
    rad_cor = shelve.open('rdd')
    vdw = shelve.open('vdw')

    df.close()
    rad_cor.close()
    vdw.close()
