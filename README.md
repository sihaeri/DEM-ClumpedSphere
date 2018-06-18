# DEM-ClumpedSphere
An algorithm to generate clumped-sphere approximation of non-spherical particles from particle STL files (X-ray micortomography, Selective Electron Microscopy etc), to be use in Discrete Element Simulations. 

If any publication is resulted from using this code, please cite: 

 1. Optimisation of blade type spreaders for powder bed preparation in additive manufacturing using DEM simulations
    ***Haeri Sina***
    *Powder Technology Vol 321, pp. 94-104, (2017)
    http://dx.doi.org/10.1016/j.powtec.2017.08.011

**Installation**
No installation is required, however, the code relies on the following two libraries
 1. Ray/Triangle Intersection version 1.0 (https://uk.mathworks.com/matlabcentral/profile/authors/1927554-jesus-p-mena-chalco)
 2. stlTools version 1.1 (https://uk.mathworks.com/matlabcentral/profile/authors/402762-pau-mico)

**Features**
 1. Can write the clumped-sphere assembly in TECPLOT format for visualisation.
 2. Outputs the moment of inertia, volume of the assembly and the errors by comparing to direct MC calculations (this is time consuming).
 3.  Writes Lammps templates of the assembly.
 4. Possibility to control surface features, (if using large smoothing a rescaling of the sphere diameters may be required, see the paper for a discussion.
