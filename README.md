# mdlibs
C header files which help read molecular dynamics trajectories

psf.h and pdb.h are for reading PSF and PDB files, respectively. Neither has any
capability to write PSF and PDB files.

psfpdb.h facilitates reading either PSF or PDB file types, and returns more
general atom information (segment name, reside name and ID, atom name, and
charge).

dcd.h is incomplete, but primarily reads DCD files. It has limited capability to
write updated coordinates to an existing DCD file.
