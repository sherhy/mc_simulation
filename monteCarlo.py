#!/usr/bin/env python3                                                          
import numpy as np
import matplotlib.pyplot as plt
from random import random, seed
from math import sqrt, log, exp

# initial conditions
seed('mc')

class Vector:
    def __init__(self, _x, _y):
        self.x = _x
        self.y = _y

    def add(self, v):
        self.x += v.x
        self.y += v.y
        return self

    def sub(self, v):
        self.x -= v.x
        self.y -= v.y

    def __str__(self):
        return f"x: {self.x} y:{self.y}"

    def distanceTo(self, p2):
        # return sqrt((self.y - p2.y)**2 + (self.x - p2.x)**2)
        return Vector(self.x - p2.x, self.y - p2.y)

    def getMagnitude(self):
        return sqrt(self.x**2 + self.y**2)

class Particle:
    border = 0

    def __init__(self, _x, _y, _name=""):
        self.pos = Vector(_x, _y)
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

def getLJP(p1, p2):
    r = p1.pos.distanceTo(p2.pos).getMagnitude()
    if r == 0 or r > 8: return 0
    return 4*((1 / r)**12 - (1 / r)**6)

def plot(plist):
    x = [p.pos.x for p in plist]
    y = [p.pos.y for p in plist]
    plt.scatter(x, y, alpha=0.7)

def rand(n=1): 
    return (random()-.5)*n

def randomVector(factor=1): 
    return Vector(rand(factor), rand(factor))

def monteCarloSim(n):
    dlessTemp = 2.74
    Particle.border = n + 1
    particles = [Particle(i,j,f"{i}{j}")
        for i in range(-n, n+1) for j in range(-n, n+1)]

    # eight ghost cells
    ghostCells = []
    ghostCells += [Particle(p.pos.x - 2*(n+1), p.pos.y - 2*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x - 0*(n+1), p.pos.y - 2*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x + 2*(n+1), p.pos.y - 2*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x - 2*(n+1), p.pos.y - 0*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x + 2*(n+1), p.pos.y + 0*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x - 2*(n+1), p.pos.y + 2*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x + 0*(n+1), p.pos.y + 2*(n+1)) for p in particles]
    ghostCells += [Particle(p.pos.x + 2*(n+1), p.pos.y + 2*(n+1)) for p in particles]
    
    N = len(particles)

    # plot(particles)

    LJP = sum([getLJP(p1, p2) for p1 in particles for p2 in particles])
    LJP += sum([getLJP(p, gp) for gp in ghostCells for p in particles])

    nudgeCount = 0
    monteCarloCycle = 0
    oldLJP = 1
    while True:
        # logic for breaking out of cycle
        if monteCarloCycle > 20000: break
        if monteCarloCycle % 1000 == 0:
            print(f"{monteCarloCycle//1000} LJP: {LJP}")
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
            particles[randIndex].pos.y)

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
        ghostCells[0*N + randIndex].pos = newLocation.add(Vector(-2*(n+1), -2*(n+1)))
        ghostCells[1*N + randIndex].pos = newLocation.add(Vector(-0*(n+1), -2*(n+1)))
        ghostCells[2*N + randIndex].pos = newLocation.add(Vector(+2*(n+1), -2*(n+1)))
        ghostCells[3*N + randIndex].pos = newLocation.add(Vector(-2*(n+1), -0*(n+1)))
        ghostCells[4*N + randIndex].pos = newLocation.add(Vector(+2*(n+1), +0*(n+1)))
        ghostCells[5*N + randIndex].pos = newLocation.add(Vector(-2*(n+1), +2*(n+1)))
        ghostCells[6*N + randIndex].pos = newLocation.add(Vector(+0*(n+1), +2*(n+1)))
        ghostCells[7*N + randIndex].pos = newLocation.add(Vector(+2*(n+1), +2*(n+1)))

        # check that energy state is conserved
        LJP += (newPotential - oldPotential)

    print(f"nudge ratio: \
        {int(nudgeCount/monteCarloCycle * 100)}% out of {monteCarloCycle}")
    return particles
    

one = monteCarloSim(5)

def g(plist, r):
    dr = 0.2
    for p in plist:
        distances = list(map(lambda point: 
            point.pos.distanceTo(p.pos).getMagnitude(), plist))
        inRange = list(filter(lambda length: length - r < dr, distances))
        inRange.remove(0) #remove the count for itself


        print(inRange)
        break

g(one, 2)
plot(one)
plt.show()

##################

# p1 = Particle(1, 1)
# p2 = Particle(0, 0)
# print(p1.pos.distanceTo(p2.pos).getMagnitude())
# print(getLJP(p1, p2))