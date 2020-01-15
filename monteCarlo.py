#!/usr/bin/env python
import shelve
from random import random, seed
from vector import Vector, randomVector
from particle import Particle

# initial conditions
seed('mc')

def getLJP(p1, p2):
    r = p1.pos.distanceTo(p2.pos).getMagnitude()
    if r == 0 or r > 8: return 0
    return 4*((1 / r)**12 - (1 / r)**6)

def makeGhost(x, y, z, plist):
    boundary = 2*(Particle.border)
    return [Particle(
        p.pos.x + x * boundary, 
        p.pos.y + y * boundary, 
        p.pos.z + z * boundary) for p in plist]

def monteCarloSim(n=2, shelved=False, stopAt=10000):
    if not shelved:
        dlessTemp = 2.74
        Particle.border = n + 1
        nudgeCount = 0
        monteCarloCycle = 0

        particles = [Particle(i,j,k,f"{i}{j}") 
            for i in range(-n, n+1) 
            for j in range(-n, n+1) 
            for k in range(-n, n+1)]

        ghostCells = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if i == 0 and j == 0 and k == 0: continue
                    ghostCells += makeGhost(k, j, i, particles)

        LJP = sum([getLJP(p1, p2) for p1 in particles for p2 in particles])
        LJP += sum([getLJP(p, gp) for gp in ghostCells for p in particles])
    else:
        n = shelved["n"]
        dlessTemp  = shelved["dlessTemp"]
        nudgeCount = shelved["nudgeCount"]
        monteCarloCycle = shelved["cycleCount"]
        particles  = shelved["plist"]
        ghostCells = shelved["ghost"]
        LJP = shelved["LJP"]

    Particle.border = n + 1
    N = len(particles)
    oldLJP = 1
    while True:
        # if nudgeCount == 1: break
        if monteCarloCycle % 1000 == 0:
            if monteCarloCycle > stopAt: break
            print(f"{monteCarloCycle//1000} LJP: {LJP}")
            db[f"{monteCarloCycle//1000}"] = {
                "n": n,
                "dimension": 3,
                "dlessTemp": dlessTemp,
                "cycleCount": monteCarloCycle,
                "nudgeCount": nudgeCount,
                "LJP": LJP,
                "plist": particles,
                "ghost": ghostCells
            }
            if abs((oldLJP - LJP)/oldLJP) < 1e-3: break
            oldLJP = LJP

        monteCarloCycle += 1

        # alter position of single atom
        randIndex = int(random()*N)
        oldPotential = sum([getLJP(particles[randIndex], p) for p in particles])
        oldPotential += sum([getLJP(particles[randIndex], p) for p in ghostCells])

        nudge = randomVector()
        newParticle = Particle(
            particles[randIndex].pos.x, 
            particles[randIndex].pos.y,
            particles[randIndex].pos.z)

        newParticle.pos.add(nudge)
        newParticle.checkLimits()

        # get potential
        particlesCopy = particles.copy()
        particlesCopy.pop(randIndex)
        newPotential = sum([getLJP(newParticle, p) for p in particlesCopy])
        newPotential += sum([getLJP(newParticle, p) for p in ghostCells])

        # print(nudge, oldPotential, newPotential)

        # move individual atoms according to LJ potential
        if newPotential > oldPotential:
            probability = exp(-(1/dlessTemp) * (newPotential - oldPotential))
            # print(f"{probability} {(newPotential - oldPotential)}")
            if random() > probability: continue
        nudgeCount += 1
        particles[randIndex].pos.add(nudge)
        particles[randIndex].checkLimits()
        newLocation = particles[randIndex].pos

        adjustedVectors = [Vector(k*2*(n+1), j*2*(n+1), i*2*(n+1))
            for k in range(-1, 2)
            for j in range(-1, 2)
            for i in range(-1, 2)
        ]
        adjustedVectors.pop(13)

        for i in range(26):
            ghostCells[i*N].pos = newLocation.add(adjustedVectors[i])

        # check that energy state is conserved
        LJP += (newPotential - oldPotential)

    print(f"nudge ratio: \
        {int(nudgeCount/monteCarloCycle * 100)}% out of {monteCarloCycle}")
    return particles
    

if __name__ == '__main__':
    db = shelve.open("mcSimulation")
    
    monteCarloSim(shelved=db['30'], stopAt=30001)

    db.close()