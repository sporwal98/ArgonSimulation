"""
CMod Project B: Periodic Boundary Conditions
This program contains methods to deal with simulation of periodic boundary conditions and images

Author: S. Porwal s1705173, F. McDougall s1713262
Version: 10/2019
"""
import numpy as np

def imageInCube(x, l):
	'''
	Image of x in the cube defined by 0<=x,y,z<=l
	
	:param x: (numpy array) vector (x_1, x_2, x_3)
	:return: image of x in cube 0<=x,y,z<=l
	'''
	return np.mod(x,l)

def closestToOrigin(x,l):
	'''
	Image of x closest to the origin, with boundary defined by cube 0<=x,y,z<=1

	:param x: (numpy array) vector (x_1, x_2, x_3)
	:return: image of x closest to the origin
	'''
	close = imageInCube(x,l)
	close[close>(l/2)] -= l
	return close
	

