"""
 Particle3D, a class to describe 3D particles

 Author: S. Porwal s1705173, F. McDougall s1713262
"""
import numpy as np
class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    position(numpy vector) - position in space
    velocity(numpy vector) - velocity in space
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, pos, vel, mass,lab):
        """
        Initialise a Particle3D instance
        
        :param pos: position as numpy vector
        :param vel: velocity as numpy vector
        :param mass: mass as float
	:param lab: label as string
        """
        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = lab
    

    def __str__(self):
        """
        Define output format.
        For particle labelled "p" at location (0.1,0.2,0.3) this will return
	"p 0.1 0.2 0.3"
        """
        return self.label + " " + str(self.position[0]) + " " +str(self.position[1])+ " " +str(self.position[2])

    
    def kineticEnergy(self):
        """
        Return kinetic energy as
        1/2*mass*|vel|^2
        """
        return 0.5*self.mass*np.linalg.norm(self.velocity)**2
        

    # Time integration methods
    def leapVelocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as numpy vector
        """
        self.velocity += dt*force/self.mass


    def leapPos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position += dt * self.velocity


    def leapPos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as numpy vector
        """
        self.position += dt * self.velocity + 0.5*dt**2*force/self.mass

    @staticmethod
    def particleFromFile(fhandle):
        '''
	Assuming line in file in following format:
	x-pos y-pos z-pos x-vel y-vel z-vel mass label

	:param fhandle: filehandle of file containing particle properties
	:return: Particle3D object with provided properties
	'''

	#Initialise Particle3D object
        line = fhandle.readline().split()

        if(len(line)==0):
            return ''
	
        p = Particle3D( np.array([float(line[0]),float(line[1]),float(line[2])]) , np.array([float(line[3]),float(line[4]),float(line[5])]) , float(line[6]) , line[7])

        return p

    @staticmethod
    def sep(p1,p2):
        '''
        Return relative vector separation of two particles
        
        :param p1: Particle3D object
        :param p2: Particle3D object
        :return: vector separation of p1 and p2: p2.pos - p1.pos
        '''
        return p2.position - p1.position
        








