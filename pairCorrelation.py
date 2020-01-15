import shelve
from particle import Particle
from vector import Vector

def g(plist, r):
    dr = 0.2
    for p in plist:
        distances = list(map(lambda point: 
            point.pos.distanceTo(p.pos).getMagnitude(), plist))
        inRange = list(filter(lambda length: length - r < dr, distances))
        inRange.remove(0) #remove the count for itself


        print(inRange)
        break

if __name__ == '__main__':
    db = shelve.open("mcSimulation")

    g(db['41']['plist'], 2)
    
    db.close()