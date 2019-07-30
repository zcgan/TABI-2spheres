# TABI-2spheres
This repo provides a version of the TABI code for the case of two dielectric spheres with multiple charges inside each sphere and surrounded by ionic solvent. The code solves the linear Poisson-Boltzmann equation for the surface potential using the Juffer et al. (1991) integral formulation and a boundary element method. The spheres are discretized using the icosahedral triangulation. The code is written in Fortran. 

# how to run it
Type make, then run the executable file "./tabipb.exe".

# input files
The parameters for TABI are in usrdata.in and the data for the two spheres is in the folder 2sphere_data.

# usrdata.in (parameters for TABI)
epsp: dielectric constant inside spheres

epsw: dielectric constant of solvent

bulk_strength: ionic concentration strength (M)

den: integer number of refinement levels on each sphere; for good accuracy chose den at least 3

order: treecode multipole expansion order; default value is order=3

maxparnode: maximum particles per leaf for octree; default value is maxparnode=500

mac: multipole acceptance criterion for treecode; default value is mac=0.5

# 2sphere_data (data for 2 spheres)
This folder contains two files, oneb_sph1.pqr, oneb_sph2.pqr, with the data describing the two spheres. Each line in the files corresponds to a charge. There are 5 columns for each line containing information about the charge as follows. All lengths are in units of angstroms and charges are in units of e.

column 1: x-coordinate of charge

column 2: y-coordinate of charge

column 3: z-coordinate of charge

column 4: charge value of charge

column 5: radius of sphere enclosing the charge

We assume each sphere has a central charge, which is the first line of each pqr file. Starting from the second line, we list the off-center charges. If there is no central charge, then set the charge value of the first line to zero.

# output file 1: energy.dat
This gives the total electrostatic free energy (Coulomb and solvation) of the 2-sphere system in units of kcal/mol. 

To find the interaction energy, compute and subtract from the above value the total electrostatic free energy for the 2-sphere system when the spheres are far from each other, say a few Debye lengths. This just need to be calculated once.

# output file 2: surface_potential.dat
This gives the surface potenial values in a format suitable for plotting with VMD.

The 1st line gives the number of vertices and the number of triangles for the entire 2-sphere system.

The next lines gives infomation about each vertex (index of vertex, xyz coordinates of vertex,	normal	vector at vertex,	surface	potential at vertex	[kcal/mol/ec], surface potential normal	derivative at vertex [kcal/mol/ec/Å].	

The next lines give the connectivity data	for	the surface	triangulation.	

# Reference
W.H. Geng	 and	 R.	 Krasny,	 A	treecode-accelerated	 boundary	 integral	 Poisson-Boltzmann	 solver	 for	
electrostatics	of	solvated	biomolecules,	J.	Comput.	Phys.,	247,	62-78	(2013).

This work was funded by National Science Foundation grants DMS-1418966 and DMS-1819094.
