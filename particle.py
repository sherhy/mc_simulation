from vector import Vector


class Particle:
    border = 3

    def __init__(self, _x, _y, _z):
        self.pos = Vector(_x, _y, _z)

    def check_limits(self):
        if abs(self.pos.x) > Particle.border or abs(self.pos.y) > Particle.border or abs(self.pos.z) > Particle.border:
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
