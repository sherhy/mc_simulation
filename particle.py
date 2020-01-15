from vector import Vector


class Particle:
    border = 0
    n = 0

    def __init__(self, _x, _y, _z, _name=""):
        self.pos = Vector(_x, _y, _z)
        self.name = _name

    def check_limits(self):
        if self.pos.x < -Particle.border:
            self.pos.x += Particle.border * 2
        elif self.pos.x > Particle.border:
            self.pos.x -= Particle.border * 2

        if self.pos.y < -Particle.border:
            self.pos.y += Particle.border * 2
        elif self.pos.y > Particle.border:
            self.pos.y -= Particle.border * 2

        if self.pos.z < -Particle.border:
            self.pos.z += Particle.border * 2
        elif self.pos.z > Particle.border:
            self.pos.z -= Particle.border * 2
