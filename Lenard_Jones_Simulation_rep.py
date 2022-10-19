#!/usr/bin/env python3
# -*- coding: utf-8 -*
"""
Created on Tue Feb  1 13:59:14 2022

"""

import numpy as np
from particle3D_copy import Particle3D as p3d
from pbc import pbc
from pbc import mic
from math import pi
import datetime
import mdutilities_20210131
import sys
import matplotlib.pyplot as plt



def force_LJ(r_ij, mod_r_ij, l, n_particles):
    '''
    Parameters
    ----------
    r_ij : Sepearation between all particles (32, 32, 3)

    sigma : defines the length constant of the lenard jones potential

    Returns
    -------
    sum_force : The net force acting on a particle as a 3D vector

    '''
    n = n_particles
    forces = np.zeros((n, n, 3))
    # sum_force = np.zeros((3,n))
    # setting the separations between all particles using the pair_sep function
    # r_ij = pair_sep(p_ij, l, n)
    for i in range(n):
        for j in range(i+1,n):
            # computing the force of the upper triangle of the matrix and then
            # just copying the negative of that on the lower triangle of the matrix
            
            # using MIC (minimum image convention) to find the nearest neighbours
            #mod_r_ij = np.linalg.norm(r_ij[i,j])
            if mod_r_ij[i,j] < 3.5:
                forces[i,j] = 48*(mod_r_ij[i,j]**-14
                                  -0.5*mod_r_ij[i,j]**-8)*(r_ij[i,j])
            else:
                forces[i,j] = 0
            forces[j,i] = - forces[i,j]
    # calculating the net force as 3D vector on every particle
    sum_force = np.sum(forces, axis= 1)
    return sum_force

def pair_sep(p_ij, l,n_particles):
    '''
    Parameters
    ----------
    p_ij : list of particle 3d instance
    
    l : side length of the lattice

    Returns
    -------
    seps : (32, 32, 3) array of the separations between particles using MIC
    '''  
    n = n_particles
    seps = np.zeros((n, n, 3))
    pos = np.zeros((n,3))
    for i, p in enumerate(p_ij):
        pos[i] = p.pos
    
    rij = pos.reshape(n,1,3) - pos.reshape(1,n,3)
    seps = mic(rij, l)
    
    '''
    for i in range(n):
            for j in range(i+1,n):
                # computing the force of the upper triangle of the matrix and then
                # just copying the negative of that on the lower triangle of the matrix

                r_ij = p_ij[i].pos - p_ij[j].pos
                # using MIC to find nearest neighbours
                seps[i,j] = mic(r_ij,l)
                seps[j,i] = - seps[i,j]
                '''
    mods = np.linalg.norm(seps, axis=2)

    return seps, mods


def pot_energy(r_ij ,l ,n_particles):
    '''
    Using symmetry in numpy arrays, this function returns the total potential energy of the system
    Parameters\
    ----------
    r_ij : Sepearation between all particles (32, 32, 3)
    
    l : side length of unit cell
    
    n_particles : total number of particles in system

    Returns
    -------
    pot_energies: total potential energy in reduced units
    '''
    n = n_particles
    pot_energy = 0.0
    # r_ij = pair_sep(p_ij, l,n)
    for i in range(n):
        for j in range(i+1, n):
            # computing the force of the upper triangle of the matrix and then
            # just copying the negative of that on the lower triangle of the matrix
            # modulus of separation
            # mod_r_ij = np.linalg.norm(r_ij[i,j])
            
            a = 1/r_ij[i,j]
            if r_ij[i,j] < 3.5:
                pot_energy += 4*(a**12-a**6)
            else:
                pot_energy += 4*((1/3.5)**12-(1/3.5)**6)
    # total_pot_energy = np.sum(np.array(pot_energies))
    
    return pot_energy


def MSD_func(Position_list, numstep, n_particles,l):
    '''
    Parameters
    ----------
    Position_list : Position of each particle as a 3d vector at each timestep. 
                    Format (numstep+1,n_particles, 3)
                    
    numstep : number of time the intergration is performed.
    
    n_particles : total number of particles in system

    Returns
    -------
    MSD : Numpy array of the mean square displacement of the particles in the simulation
    '''

    p_0 = Position_list[0]
    MSD = np.zeros(numstep+1)
    for i in range(numstep+1):
        MSD[i] = (np.linalg.norm(mic(Position_list[i]-p_0, l)))**2
    
    MSD = MSD/n_particles
    return MSD



def RDF_func(Seps_list, l, n_particles, rho):
    '''
    Parameters
    ----------
    Seps_list : Separations between of each particle as a 3d vector at each timestep.
                Format (numstep+1,n_particles, n_particles, 3)
    
    l : side length of unit cell
    
    rho : density of system
    
    n_particles : total number of particles in system

    Returns
    -------
    g(r) : y axis of the RDF histogram
    
    r : x axis of the RDF histogram
    ''' 
    # Setting up fixed bins for the histogram, so that I can compare between different runs
    # sqrt(3)/2 l is because that is the maximum
    Bins = np.linspace(0, l * ((3)**0.5)/2, 100)
    # making the separation vectors as magnitudes
    # mag_seps = np.linalg.norm(Seps_list,axis=-1)
    mag_seps = Seps_list

    flat_mag_seps=mag_seps.flatten()

    # creating my unweighted g(r) and the bin edges as,
    # ignoring the sep = 0, as those are not imprtant in this case
    g_R, R_list = np.histogram(flat_mag_seps[flat_mag_seps>0], bins = Bins)
    
    # bin size
    dr = R_list[1] - R_list[0]
    
    # using the middle of bins to find the midpoint between bins to use in plotting
    r = R_list[1:]-0.5*dr
    
    # creating the weights for the RDF function
    N_p0 = (n_particles * 4*pi * r**2 * rho * dr)*Seps_list.shape[0]

    
    return g_R/N_p0, r



def main():
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        sys.exit(-1)
    else:
        outfile_name = sys.argv[1]
        
    outfile = open(outfile_name, "w")
    Vel_file = open('Velocities_file', 'w')
    Energies_file = open('Energies_file', 'w')
    Energies_file.write('Time [A(amu/eV)^0.5], Total Energy [eV], Kinetic Energy [eV], Potential Energy [eV]\n')
    MSD_file = open('MSD_file', 'w')
    MSD_file.write('Time [A(amu/eV)^0.5], MSD\ [A^2]n')
    RDF_file = open('RDF_file', 'w')
    RDF_file.write('r [A], g(r)\n')
    
    file ='solid_argon_report.dat' #str(input('please input the name of the file to use for the run (in the solid_argon.dat format):'))
    with open(file, 'r') as f:
        f.readline()
        dt_numstep_data = f.readline().split()
        f.readline()
        sim_parameters = f.readline().split()

    #Simulation inputed parameters set up from the imported file
    k=1
    dt = float(dt_numstep_data[0])*k
    numstep= int(int(dt_numstep_data[1])/1)
    n_particles = int(sim_parameters[0])
    rho =(((float(sim_parameters[1]))/(10**6*6.0221409e23))*(3.405e-10)**-3)**(-1)
   
    print('rho =', rho)
    mytemp = 0.6 #float(sim_parameters[2])/119.8
    print('ε/kb =', mytemp)
    print('dt =', dt)
    # creating my particle3D instances
    p_ij = []
    for i in range(n_particles):
        p_ij.append(p3d(f'p{i}', 1.0 ,np.zeros(3), np.zeros(3)))
        
    # setting the initial positions of the particles
    box_length, _ = mdutilities_20210131.set_initial_positions(rho, p_ij)
    
    l = box_length[0]
    
    # setting up the initial velocities of the particles
    mdutilities_20210131.set_initial_velocities(mytemp, p_ij, seed = 123456789)
    


    #initial conditions for time intergrator set up
    p_0 = np.zeros((n_particles,3))
    for i in range(n_particles):
        p_0[i] = p_ij[i].pos

    r_ij_i, mod_rij = pair_sep(p_ij, l, n_particles)
    
    forces = force_LJ(r_ij_i, mod_rij,  l, n_particles)


    total_energy = pot_energy(mod_rij, l, n_particles) + p3d.sys_kinetic(p_ij)
  

    
    time = 0.0
    Time_list = np.zeros(numstep+1)
    Time_list[0] = time

 
    Tot_Energy_list = np.zeros(numstep+1)
    Tot_Energy_list[0]=total_energy

    Energies_file.write('{0:f}, {1:12.8f}, {2:12.8f}, {3:12.8f}\n'.format(
        Time_list[0], total_energy, p3d.sys_kinetic(p_ij), pot_energy(mod_rij, l, n_particles)))
    
    Position_list = np.zeros((numstep+1,n_particles, 3))
    Position_list[0] = p_0
    
    Seps_list = np.zeros((numstep+1,n_particles, n_particles))
    Seps_list[0] = mod_rij 

    #wrtting to file the particle initial conditions used for simulation
    outfile.write(f"{n_particles}\n")
    outfile.write(f'timestep = {0}\n')
    Vel_file.write(f"{n_particles}\n")
    Vel_file.write(f'timestep = {0}\n')
    for p in p_ij:
        outfile.write(str(p))
        Vel_file.write(p.__strvel__())
#Start of time intergration
    for i in range(numstep):
        #print(i)
        p_t = np.zeros((n_particles,3))

        for j in range(n_particles):
            
            #updating to the 2nd order the positions
            p_ij[j].update_pos_2nd(dt,forces[j])
            p_ij[j].pos = pbc(p_ij[j].pos,l)
            #applying pbc to the updated positions 
 

            p_t[j] = p_ij[j].pos


        Position_list[i+1] = p_t
        r_ij, mod_rij = pair_sep(p_ij, l, n_particles)

        forces_new = force_LJ(r_ij, mod_rij, l, n_particles)
        # update velocities using mean force between the forces calculated at t-dt and t+dt
        for j in range(n_particles):

            p_ij[j].update_vel(dt, 0.5*(forces[j]+forces_new[j]))


        # Re-define force value     
        forces = forces_new



        Seps_list[i+1] = mod_rij
        # Increase time
        time += dt
        #print(time)
        pe =  pot_energy(mod_rij, l, n_particles)
        total_energy = pe + p3d.sys_kinetic(p_ij)


        Time_list[i+1] = time

        Tot_Energy_list[i+1] = total_energy

        Energies_file.write('{0:f}, {1:12.8f}, {2:12.8f}, {3:12.8f}\n'.format(
            time, total_energy, p3d.sys_kinetic(p_ij), pe))

        #setting up the output file format
        outfile.write(f"{n_particles}\n")
        outfile.write(f'timestep = {i}\n')
        Vel_file.write(f"{n_particles}\n")
        Vel_file.write(f'timestep = {i}\n')
        # print(f'{i}, 7')
        for p in p_ij:
            outfile.write(str(p))
            Vel_file.write(p.__strvel__())

    


    MSD = MSD_func(Position_list, numstep, n_particles,l)
    
    for i in range(numstep+1):
        MSD_file.write('{0:f}, {1:12.8f}\n'.format(Time_list[i], MSD[i]))

    
    g_r, r, = RDF_func(Seps_list, l, n_particles, rho)
    
    for i in range(len(r)):
        RDF_file.write(f'{r[i]}, {g_r[i]}\n') 
      
    # Energy vs time        
    plt.title(f'E vs time {file} {dt}')
    plt.xlabel('Time [A(amu/eV)^0.5]')
    plt.ylabel('Total Energy [eV]')
    plt.plot(Time_list, Tot_Energy_list)
    plt.show()
    plt.clf()
    # MSD vs Time
    plt.title(f'MSD vs time {file} {dt}')
    plt.xlabel('Time [A(amu/eV)^0.5]')
    plt.ylabel('Mean Square Displacement [A^2]')
    plt.plot(Time_list, MSD)
    plt.show()
    plt.clf()
    # RDF vs r
    plt.title(f'g(r) vs r {file} {dt}')
    plt.xlabel('r [A]')
    plt.ylabel('g(r)')
    plt.plot(r,g_r)
    plt.show()
    

    outfile.close()
    Energies_file.close()
    MSD_file.close()
    RDF_file.close()


a=datetime.datetime.now()
if __name__ == "__main__":
    main()
b=datetime.datetime.now()
print('Run Time', b-a)
