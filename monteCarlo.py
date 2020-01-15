#!/usr/bin/env python
import shelve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from random import random, seed
from math import sqrt, log, exp

# initial conditions
seed('mc')

class Vector:
    def __init__(self, _x, _y, _z):
        self.x = _x
        self.y = _y
        self.z = _z

    def add(self, v):
        self.x += v.x
        self.y += v.y
        self.z += v.z
        return self

    def __str__(self):
        return f"x: {self.x} y:{self.y}"

    def distanceTo(self, p2):
        # return sqrt((self.y - p2.y)**2 + (self.x - p2.x)**2)
        return Vector(self.x - p2.x, self.y - p2.y, self.z - p2.z)

    def getMagnitude(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)

class Particle:
    border = 0

    def __init__(self, _x, _y, _z, _name=""):
        self.pos = Vector(_x, _y, _z)
        self.name = _name

    def checkLimits(self):
        if self.pos.x  < -Particle.border: 
            self.pos.x += Particle.border * 2
        elif self.pos.x > Particle.border: 
            self.pos.x -= Particle.border * 2

        if self.pos.y  < -Particle.border: 
            self.pos.y += Particle.border * 2
        elif self.pos.y > Particle.border: 
            self.pos.y -= Particle.border * 2

        if self.pos.z  < -Particle.border: 
            self.pos.z += Particle.border * 2
        elif self.pos.z > Particle.border: 
            self.pos.z -= Particle.border * 2

def getLJP(p1, p2):
    r = p1.pos.distanceTo(p2.pos).getMagnitude()
    if r == 0 or r > 8: return 0
    return 4*((1 / r)**12 - (1 / r)**6)

def rand(n=1): 
    return (random()-.5)*n

def randomVector(factor=1): 
    return Vector(rand(factor), rand(factor), rand(factor))

def makeGhost(x, y, z, plist):
    boundary = 2*(Particle.border)
    return [Particle(
        p.pos.x + x * boundary, 
        p.pos.y + y * boundary, 
        p.pos.z + z * boundary) for p in plist]


def monteCarloSim(n, shelved = False):
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
        if monteCarloCycle > 20000: break
        if monteCarloCycle % 1000 == 0:
            print(f"{monteCarloCycle//1000} LJP: {LJP}")
            archive[f"{monteCarloCycle//1000}"] = {
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
    
def g(plist, r):
    dr = 0.2
    for p in plist:
        distances = list(map(lambda point: 
            point.pos.distanceTo(p.pos).getMagnitude(), plist))
        inRange = list(filter(lambda length: length - r < dr, distances))
        inRange.remove(0) #remove the count for itself


        print(inRange)
        break

def plot(plist):
    x = [p.pos.x for p in plist]
    y = [p.pos.y for p in plist]
    z = [p.pos.z for p in plist]
    ax.scatter(x, y, z, alpha=0.7)
    return ax

def update(i):
    return plot(archive[f"{i}"]['plist'])


if __name__ == '__main__':
    archive = shelve.open("mcSimulation")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # one = monteCarloSim(2, archive['10'])
    # g(one, 2)
    plot(archive['0']['plist'])
    # plot(archive['1']['plist'])
    # plot(archive['2']['plist'])
    # plot(archive['3']['plist'])
    plot(archive['4']['plist'])
    # plot(archive['5']['plist'])
    # plot(archive['6']['plist'])
    # plot(archive['7']['plist'])
    plot(archive['8']['plist'])
    # plot(archive['9']['plist'])
    # plot(archive['10']['plist'])
    plot(archive['12']['plist'])
    # plot(archive['20']['plist'])
    plt.show()
    archive.close()

    #tricky to model 3d of points, just because we're not used to it;; maybe use
    #p5js or something for some rotation to see the points from rotated perpsecitves