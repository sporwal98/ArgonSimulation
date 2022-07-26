#Contains methods to operate on lists of particles and calculate Lennard-Jones Forces and Energies
import numpy as np
import Particle3D as p3d
import LennardJonesPotential as ljp


def multiForce(particles,epsilon,sigma,cutoff, box):
    #Returns a list of net force on each particle in the box
    forces = []
    for i in range(len(particles)):
        forces.append(0)
    for i in range(len(particles)):
        for j in range(i+1,len(particles)):
            force = ljp.forceLJP(particles[i].position,particles[j].position,cutoff,box)
            forces[i]+=force
            forces[j]-=force
        
    return forces

def multiPE(particles,epsilon,sigma,cutoff,box):
    #Returns the total potential energy of the system of particles
    PE = 0
    for i in range(len(particles)):
        for j in range(i+1,len(particles)):
            PE+=ljp.potEnergyLJP(particles[i].position,particles[j].position,cutoff,box)
    return PE

def multiKE(particles):
    #Returns a list of kinetic energies of each particle in the box
    KE = []
    for p in particles:
        KE.append(p.kineticEnergy())
    return np.array(KE)















