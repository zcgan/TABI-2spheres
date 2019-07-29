# TABI-2spheres
TABI specifically for the electrostatic interaction between two dielectric spheres, with multiple charges inside each sphere.

# How to run it
Simply type make, then run the executable	file by "./tabipb.exe"  (using either gfortran or ifort compiler)

# about input files
# inputfile #1: usrdata.in
epsp: the dielectric constant inside the spheres
epsw: the dielectric constant of the solvent
bulk_strength: the ionic concentration (M)
den: an integer that determines the number of grids on each sphere (20*4^den)
order: treecode multipole expansion order
maxparnode: the maximum particles per leaf for the tree structure
mac: multipole acceptance criterion for treecode

# inputfile #2: \2sphere_data\oneb_sph1.pqr and \2sphere_data\oneb_sph2.pqr
