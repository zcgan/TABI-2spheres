# TABI-2spheres
TABI specifically for the electrostatic interaction between two dielectric spheres, with multiple charges inside each sphere.

# How to run it
Simply type make, then run the executable	file by "./tabipb.exe"  (using either gfortran or ifort compiler)

# about input files
# inputfile #1: usrdata.in
epsp: the dielectric constant inside the spheres

epsw: the dielectric constant of the solvent

bulk_strength: the ionic concentration (M)

den: an integer that determines the number of icosahedron grids on each sphere (20*4^den), for a good accuracy, choosing den at least 3.

order: treecode multipole expansion order, also need to be an integer, here set to be 3.

maxparnode: the maximum particles per leaf for the tree structure

mac: multipole acceptance criterion for treecode

# inputfile #2: \2sphere_data\oneb_sph1(2).pqr
There are 5 columns for each pqr file for sphere #1 and #2, that contains the information for their inside charge locations and magnitudes, and radius of each sphere.

All length are in unit of angstrom, charges in unit of e.

col#1: x-coordinate for each point charge

col#2: y-coordinate for each point charge

col#3: z-coordinate for each point charge

col#4: the charge value for each point charge

col#5: the radius of the sphere that enclosing the charge

Note that for simplicity, we assume each sphere has a central charge, which is the first row of each pqr file. Starting from the second row, we list the off-centered charges. (If there is no central charge, just set the charge value of the first row to be 0)

# Output file #1: energy.dat
it gives the electrostatic free energy value in unit of Kcal/mol. 

Note that if one wants to get the interaction energy, then further subtract the energy value for the two spheres at sufficient large distance.

# Output file #2: surface_potential.dat
containing	

(1)	number	of	nodes,	number	of	triangles,	

(2)	node	index,	vertices,	normal	vectors,	surface	potentials	[kcal/mol/ ec],	surface	
potential	normal	derivatives	[kcal/mol/ec/Å],	

(3)	connectivity	data	for	icosahedron surface	triangulation.	

# Reference
 W.H. Geng	 and	 R.	 Krasny,	 A	treecode-accelerated	 boundary	 integral	 Poisson-Boltzmann	 solver	 for	
electrostatics	of	solvated	biomolecules,	J.	Comput.	Phys.,	247,	62-78	(2013).
