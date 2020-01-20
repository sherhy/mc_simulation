import math
import shelve


def calculate_pressure(mc, rdd, name):
    def ljp_prime(r, dr):
        lr = dr + r
        return 4 * (((1 / lr) ** 12 - (1 / lr) ** 6) - ((1 / r) ** 12 - (1 / r) ** 6))

    rho = mc[name]['n'] / (int(mc[name]['border']) * 2) ** 3
    kt = mc[name]['dlessTemp']
    sigma = 1
    r_cutoff = 5

    gs = rdd[name + "_g"]
    rs = rdd[name + "_r"]
    # increments = [rs[i] - rs[i - 1] for i in range(1, len(rs))]

    integral = sum([rs[i] ** 3 * ljp_prime(rs[i], rs[i] - rs[i - 1]) * gs[i] for i in range(0, len(rs) - 1)])

    pressure = rho * kt - math.pi * (2 / 3) * rho ** 2 * integral
    # constant_coefficient = (16 / 9) * math.pi * rho ** 2 * sigma ** 3
    # analytic_pressure = constant_coefficient * (2 * (sigma / r_cutoff) ** 9 - 3 * (sigma / r_cutoff) ** 3)
    # pressure += analytic_pressure

    return pressure


if __name__ == "__main__":
    reduced_volume = 0.75
    db_name = '0.75_7'
    with shelve.open("./db/mc") as db:
        with shelve.open("./db/rdd") as rad:
            with shelve.open("./db/vdw") as vdw:
                vdw[db_name] = [reduced_volume, calculate_pressure(db, rad, db_name)]
