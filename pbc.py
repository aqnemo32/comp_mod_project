# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 15:06:44 2021

@author: achil
"""


import numpy as np

def pbc(r_12, l):
    '''

    Parameters
    ----------
    x : 3 dimensional vector in the form of a numpy array
        
    l : side length of cubc unit cell with one corner as the 0, 0, 0 point


    Returns
    -------
    the image of x in the unitcell of side length l and with a corner as the 0, 0, 0 point


    '''
    #taking the modulo of the original position vector wrt to l
    y=np.mod(r_12 , l)

    return y

def mic(r_12, l):
    '''
    
    Parameters
    ----------
    x : 3 dimensional vector in the form of a numpy array

    l : side length of cubc unit cell with one corner as the 0, 0, 0 point


    Returns
    -------
    The image of x closest to the origin by using the np.where() function


    '''
    #creating the image of x in unit cell in case user input of x is not in cell of side length l
    #finding the image of x closest to the origin (0,0,0)----In order to do so I first check if the coordinates in y-->
    #--> are bigger or smaller than l/2. This is where the function np.where() is useful. The function checks wether-->
    #--> a coordinate in y is bigger than l/2, if bigger the y coordinate is substracted by l if not remains the same
    #eg for l=4 and y=(1,3,3) the image closest to the origin would be y'=(1,3-l,3-l) or (1,-1,-1)
    closest_image= np.mod(r_12+l/2,l)-l/2
    return closest_image

def main():
    #user input of the original position vector
    x_0= float(input('Input value for x_0: '))
    x_1= float(input('Input value for x_1: '))
    x_2= float(input('Input value for x_2: '))
    
    #user input of the side length of unit cell
    l= float(input('Input value for l: '))
    
    # creating the x array used in the two functions by using the user inputed values
    x=np.array([x_0,x_1,x_2])
    
    pbc(x,l)
    
    mic(x,l)


if __name__ == "__main__":
    main()