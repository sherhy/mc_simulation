#!/usr/bin/env python
import shelve
from math import exp
from random import random

from particle import Particle
from vector import Vector


def rand(n=1):
    return (random() - .5) * n


def random_vector(factor=1):
    return Vector(rand(factor), rand(factor), rand(factor))


def get_ljp(p1, p2):
    r = p1.pos.get_distance_to(p2.pos)
    if r == 0 or r > 8:
        return 0
    return 4 * ((1 / r) ** 12 - (1 / r) ** 6)


def create_ghost(x, y, z, plist):
    boundary = 2 * Particle.border
    return [Particle(
        p.pos.x + x * boundary,
        p.pos.y + y * boundary,
        p.pos.z + z * boundary) for p in plist]


def run_monte_carlo(db, n=125, kt=2.74, shelved: dict = False, stop_at=10, dont_stop=False):
    if shelved:
        n = shelved["n"]
        Particle.border = shelved["border"]
        dless_temp = shelved["dlessTemp"]
        nudge_count = shelved["nudgeCount"]
        monte_carlo_cycle = shelved["cycleCount"]
        particles = shelved["plist"]
        ghost_cells = shelved["ghost"]
        total_ljp = shelved["LJP"]
    else:
        dless_temp = kt
        nudge_count = 0
        monte_carlo_cycle = 0

        particles = [Particle(i, j, k) for i in range(-2, 3) for j in range(-2, 3) for k in range(-2, 3)]
        if 125 < n:
            for _ in range(abs(n - 125)):
                particles.append(Particle(rand(2), rand(2), rand(2)))
        elif 125 > n:
            for _ in range(abs(n - 125)):
                rand_index = int(random() * len(particles))
                particles.pop(rand_index)

        ghost_cells = list()
        for k in range(-1, 2):
            for j in range(-1, 2):
                for i in range(-1, 2):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    ghost_cells += create_ghost(i, j, k, particles)

        total_ljp = sum([get_ljp(p1, p2) for p1 in particles for p2 in particles])
        total_ljp += sum([get_ljp(p, gp) for gp in ghost_cells for p in particles])

    min_run = 4
    Particle.n = len(particles)
    Particle.reduced_volume = (2 * Particle.border) ** 3 / Particle.n
    old_total_ljp = 1
    while True:
        # if nudge_count == 1: break
        if monte_carlo_cycle % 1000 == 0:
            cycle = monte_carlo_cycle // 1000
            print(f"{Particle.reduced_volume:.2f}_{cycle} LJP: {total_ljp}")
            db[f"{Particle.reduced_volume:.2f}_{cycle}"] = {
                "n": n,
                "border": Particle.border,
                "dimension": 3,
                "dlessTemp": dless_temp,
                "cycleCount": monte_carlo_cycle,
                "nudgeCount": nudge_count,
                "LJP": total_ljp,
                "plist": particles,
                "ghost": ghost_cells,
            }
            db[f"{Particle.reduced_volume:.2f}"] = cycle
            if cycle >= stop_at:
                break
            elif cycle <= min_run or dont_stop:
                pass
            elif abs((old_total_ljp - total_ljp) / old_total_ljp) < 1e-3:
                break
            old_total_ljp = total_ljp

        monte_carlo_cycle += 1

        # alter position of single atom
        rand_index = int(random() * Particle.n)
        old_potential = sum([get_ljp(particles[rand_index], p) for p in particles])
        old_potential += sum([get_ljp(particles[rand_index], gp) for gp in ghost_cells])

        nudge = random_vector()
        new_particle = Particle(particles[rand_index].pos.x, particles[rand_index].pos.y, particles[rand_index].pos.z)
        new_particle.pos.add(nudge)
        new_particle.check_limits()

        # get potential
        particles_copy = particles.copy()
        particles_copy.pop(rand_index)  # ljp explodes where r << 1
        new_potential = sum([get_ljp(new_particle, p) for p in particles_copy])
        new_potential += sum([get_ljp(new_particle, gp) for gp in ghost_cells])

        # move individual atoms according to LJ potential
        if new_potential > old_potential:
            probability = exp(-(1 / dless_temp) * (new_potential - old_potential))
            if random() > probability:
                continue

        # update ghost cells
        nudge_count += 1
        particles[rand_index].pos.add(nudge)
        particles[rand_index].check_limits()
        new_location = particles[rand_index].pos

        boundary = 2 * Particle.border
        adjusted_vectors = [Vector(i * boundary, j * boundary, k * boundary)
                            for k in range(-1, 2) for j in range(-1, 2) for i in range(-1, 2)]
        adjusted_vectors.pop(13)  # null vector

        for i in range(26):
            p = Particle(new_location.x, new_location.y, new_location.z)
            p.pos.add(adjusted_vectors[i])
            ghost_cells[rand_index + i * Particle.n] = p

        total_ljp += (new_potential - old_potential)

    print(f"nudge ratio: {int(nudge_count / monte_carlo_cycle * 100)}% out of {monte_carlo_cycle}")
    return monte_carlo_cycle // 1000


if __name__ == '__main__':
    mc = shelve.open("mc")

    # change two parameters
    # _from = 0
    # run_monte_carlo(mc, n=125, kt=2.74)
    run_monte_carlo(mc, shelved=mc['1.73_10'], stop_at=20)

    mc.close()
