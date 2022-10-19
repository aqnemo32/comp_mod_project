# Computer Modeling Project
The aims of this experiment were to find the ideal timestep using convergence for a Velocity Verlet simulation of gaseous and solid argon interacting via a Lennard-Jones potential.
## Instructions:

1) Run the code; command line: python3 Lenard_Jones_Simulation.py LJ_output.xyz
2) For the file name write: solid_argon.dat or gas_argon.dat
3) For the animation, open command prompt at file location; command line : python3 xyz_plot.py --3d LJ_output.xyz


### UNITS:

Time = [σ(m/ε)<sup>0.5</sup>] [A(amu/ε)<sup>0.5</sup>]

Energies = [ε] [ε]

MSD = [σ<sup>2</sup>] [A<sup>2</sup>]

r = [σ] [A]

Temp = [ε/k_b] [ε/(ε*K)]

g(r) = unit-less as probability distribution

### Input files format:

First line : comment with dt and numstep

Second line : values for dt and numstep respectively

Third line : comment with number of particles, density and temperature

Fourth line: values for number of particles, density and temperature

    # dt	Numstep
    5e-3	2000
    # N_Particles	Rho	Temperature(reduced units)
    30	0.05	1.0

Output files format:
These are all in column format

### Energies_file:

		Time, Total Energy, Kinetic Energy, Potential Energy
		0.0, 1.3, 0.6, 0.7

		MSD_file:
		Time, MSD
		0.0, 1.3

		RDF_file:
		r, g(r)
		0.0, 0.64

### LJ_output.xyz:

This format is repeated for every timestep

Number of particles

Timestep

Label x_pos y_pos z_pos 

eg:

  30
  timestep = 0
  p0  0.0   0.0   0.0
  p1  0.0   2.108581663254373   2.108581663254373 etc, etc

### Velocities_file:

This format is repeated for every timestep

Number of particles

Timestep

Label x_vel y_vel z_vel

eg:

  32
  timestep = 0
  p0  0.5401708158818058   0.5448161766866988   0.08695675334517412
  p1  -1.7165043075046114   -0.8060875966399225   -0.49719739227775783 etc, etc
