#Main simulation code for Project B: Liquid and gas simulations
import numpy as np
from Particle3D import Particle3D
import MultiParticleUtilities as mpu
import MDUtilities as mdu
import pbc
import sys
import matplotlib.pyplot as plt

def main():
    #Simulates an N-body Lennard-Jones simulation of Argon using Velocity Verlet Time Intergration
    if len(sys.argv)!= 5 :
        print('Wrong number of arguments')
        quit()
    else:
        #Setting up the files for input and output
        infile = sys.argv[1]
        trajfile = sys.argv[2]
        efile = sys.argv[3]
        rdffile = sys.argv[4]

    infile = open(infile,'r')

    #Scale factor for reduced units of Temperature
    Tstar = 119.8
    
    #Reading simulation parameters from input file
    line = infile.readline().split()
    epsilon,sigma,cutoff= float(line[0]),float(line[1]),float(line[2])
    av = 6.023e23
    rho,temp,dt,numstep = (1e6*av*sigma**3/float(line[3])),(float(line[4])/Tstar),float(line[5]),int(line[6])
    stepobs,m,N,nbins = int(line[7]),float(line[8]),int(line[9]), int(line[10])

    infile.close()
    
    #Opening output files
    trajfile = open(trajfile,'w')
    efile = open(efile,'w')
    efile.write('#Time\tEnergy\tKE\tPE\tMSD(t)\n')

    time = 0
    step = 0

    #Initialise list of particles
    particles = []
    for i in range(N):
        particles.append(Particle3D(np.array([0,0,0]),np.array([0,0,0]),1,('p'+str(i))))
        

    #Initialising particle positions and velocities
    box = mdu.set_initial_positions(rho, particles)
    vmdout(particles,trajfile,time)
    mdu.set_initial_velocities(temp,particles)
    forces = mpu.multiForce(particles, epsilon, sigma, cutoff,box[0])

    times = []
    msds = []
    rdf = []

    originals = []
    for p in particles:
        originals.append(p.position)


    while(step<numstep):
        #time integration
        
        time+=dt
        for i in range(len(particles)):
            #Using second order position update to update positions of particles
            particles[i].leapPos2nd(dt,forces[i])

        for p in particles:
            #Ensuring all particles are within the periodic boundary conditions
            p.position = pbc.imageInCube(p.position,box[0])
        
        #Updating forces
        newforces = mpu.multiForce(particles,epsilon,sigma,cutoff,box[0])

        for i in range(len(particles)):
            #Updating velocities
            particles[i].leapVelocity(dt, 0.5*(forces[i]+newforces[i]))

        forces = newforces
        

        if(step%stepobs == 0):
            #Writing observables
            #print(step)
            vmdout(particles,trajfile,time)
            times.append(time)
            msds.append(msd(particles,originals,box[0]))
            E = mpu.multiPE(particles,epsilon,sigma,cutoff,box[0]) + np.sum(mpu.multiKE(particles))
            wr = ('{:.4f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:f}\n').format(time,E,(np.sum(mpu.multiKE(particles))), (mpu.multiPE(particles,epsilon,sigma,cutoff,box[0])),msds[-1]) 
            efile.write(wr)
            rdf.append(rdf_calc(particles,rho,box[0]))
        

        step += 1
    
    #Creating a histogram from rdf data
    hist, bin_eds = np.histogram(rdf, bins = nbins,density = True)

    #Separating r and rdf
    bin_cents = []
    for i in range(len(bin_eds)-1):
        bin_cents.append((bin_eds[i]+bin_eds[i+1])/2)
    

    rdffile = open(rdffile,'w')
    rdffile.write('#r\tRDF(r)\n')
    #Outputting rdf to file
    for i in range(len(bin_cents)):
        rdffile.write('{:.2f}\t{:.2f}\n'.format(bin_cents[i],hist[i]))
    efile.close()
    rdffile.close()
    trajfile.close()

    #Normalising the histogram values
    for i in range(len(bin_eds)-1):
        hist[i]/=rdf_norm_fact(bin_cents[i],rho,(bin_eds[i+1]-bin_eds[i]))

    #Plotting the rdf as a function of r
    plt.figure()
    plt.title('RDF vs r')
    plt.plot(bin_cents,hist)
    plt.xlabel('r')
    plt.ylabel('rdf')
    plt.show()

    #Plotting the MSD as function of time
    plt.figure()
    plt.title('MSD vs time')
    plt.plot(times,msds)
    plt.xlabel('t')
    plt.ylabel('msd')
    plt.show()
    
    

def msd(particles,originals,boxsize):
    #Returns Mean-Squared Displacement of particles from original positions
    sums = 0
    for i in range(len(particles)):
        sums = sums + (np.linalg.norm(pbc.closestToOrigin(particles[i].position - originals[i], boxsize)))**2
    return sums/len(particles)

def rdf_calc(particles,rho,boxsize):
    #Returns list of distances between particles
    dists = []
    for i in range(len(particles)):
        for j in range(i+1,len(particles)):
            dists.append(np.linalg.norm(pbc.closestToOrigin(Particle3D.sep(particles[i],particles[j]),boxsize)))
    for i in range(len(dists)):
        dists[i]/=(len(particles)*rho)
    return dists

def rdf_norm_fact(r,rho,histdiff):
    #Returns the normalisation factor for the RDF
    return (4*np.pi*r*r*rho*histdiff)

def vmdout(particles,trajfile,time):
    #Outputs current location of all particles in VMD compatible format
    trajfile.write(str(len(particles))+'\n')
    trajfile.write('Time = '+str(time))
    for p in particles:
        trajfile.write('\n'+str(p))
    trajfile.write('\n')

if __name__ == "__main__":
    main()
