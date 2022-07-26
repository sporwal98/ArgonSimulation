import numpy as np
import pbc

def potEnergyLJP(r1,r2,cutoff,box):
    #Returns the Potential Energy of two particles at r1 and r2 interacting in Lennard-Jones Potential
    r = np.linalg.norm(pbc.closestToOrigin(r2 - r1,box))
    if(r>cutoff):
        return 0
    potE = 4*(np.power((1/r),12)-np.power((1/r),6))
    return potE

def forceLJP(r1,r2,cutoff,box):
    #Returns force on particle 1 at r1 due to particle 2 at r2 interacting with Lennard-Jones Potential
    rv = pbc.closestToOrigin(r1-r2, box)
    r = np.linalg.norm(rv)
    if(r>cutoff):
        return np.array([0.,0.,0.])
    force = 48*((1/np.power(r,14)) - (1/(2*np.power(r,8)))) * rv
    return force
