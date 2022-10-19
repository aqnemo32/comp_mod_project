"""
 CompMod Ex2: Particle3D, a class to describe point particles in 3D space

 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

 Includes time integrator methods +...


"""
import numpy as np
import sys


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_p3d - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        self.label = label
        
        self.mass = float(mass)
        self.vel = np.array(velocity, float)
        self.pos = np.array(position, float)
        
#checking if the iniput of mass is in the correct type (ie float)
#if not the exit(-1) is applied and code stop runing to prevent it crashing
        if type(mass) != type(0.1):
            print('error')
            sys.exit(-1)
#checking if the dimensions of pos and vel arrays are correct
#If they are not correct the exit(-1) is applied and the code stops running
        if len(position) != 3:
            print('error')
            sys.exit(-1)
        

        if len(velocity) != 3:
            print('error')
            sys.exit(-1)


    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        """
        xyz_string = f'{self.label}  {self.pos[0]}   {self.pos[1]}   {self.pos[2]}\n'
        return xyz_string
    
    def __strvel__(self):
        return f'{self.label}  {self.vel[0]}   {self.vel[1]}   {self.vel[2]}\n'


    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        """
        ke = 0.5 * self.mass * np.linalg.norm(self.vel)**2
        return ke


    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        
        :return momentum: np.array, m*v
        """
        momentum = self.mass * self.vel
        return momentum


    def update_pos(self, dt):
        """
        Updates the position of a particle to 1st order,
        r(t + dt) = r(t) + dt ·v(t),
        where the timestep dt is given by dt.
        """
        dt = float(dt)
        self.pos = self.pos + self.vel * dt


    def update_pos_2nd(self, dt, force):
        """
        Updates the position of a particle to 2nd order,
        r(t + dt) = r(t) + dt ·v(t) + dt2 ·f(t)/2m,
        where the timestep dt is given by dt, and the force f by force
        """
        self.dt = float(dt)
        self.force = np.array(force, float)
        self.pos = self.pos + self.dt * self.vel + (self.dt)**2 * self.force/(2 * self.mass)


    def update_vel(self, dt, force):
        """
        Updates the velocity of a particle to 1st order,
        v(t + dt) = v(t) + dt ·f(t)/m
        where the timestep dt is given by dt, and the force f by force
        """
        self.dt = float(dt)
        self.force = np.array(force, float)
        self.vel = self.vel + self.dt * self.force/self.mass


    @staticmethod
    def new_p3d(file_handle):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one line per planet in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param inputFile: Readable file handle in the above format

        :return Particle3D instance
        """
        a=file_handle.readline().split()
        label = a[0]
        mass = float(a[1])
        position = np.array([a[2],a[3],a[4]], float)
        velocity = np.array([a[5],a[6],a[7]], float)
        return Particle3D(label, mass, position, velocity)


    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Computes the sum of the kinetic energies of all the particles in p3d_list
        
        :param p3d_list: list in which each item is a P3D instance
        :return sys_ke: the total kinetic energy of the system 
        """
        sys_ke = 0
        
        for i in range(len(p3d_list)):
            sys_ke += 0.5*p3d_list[i].mass * np.linalg.norm(p3d_list[i].vel)**2
        return sys_ke


    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system of particles 
        :return com_vel: Centre-of-mass velocity, total mass of system
        """
        # Grabs the mass data from the p3d_list
        total_mass = 0
        for i in range(len(p3d_list)):
            total_mass += p3d_list[i].mass
        #center of mass velocity
        com_vel = 0
        for i in range(len(p3d_list)):
            com_vel += (p3d_list[i].momentum())/total_mass
        return total_mass, com_vel
