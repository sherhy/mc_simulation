from math import sqrt
from random import random


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

    def get_distance_to(self, p2):
        # return sqrt((self.y - p2.y)**2 + (self.x - p2.x)**2)
        return Vector(self.x - p2.x, self.y - p2.y, self.z - p2.z)

    def get_magnitude(self):
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)


def rand(n=1):
    return (random() - .5) * n


def random_vector(factor=1):
    return Vector(rand(factor), rand(factor), rand(factor))
