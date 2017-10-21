#ifndef LATTICE_H
#define LATTICE_H

struct latt //structure that contains data about direct and reciprocal lattice (same as defined in vasp, lattice.inc)
{
        double SCALE;
        double A[3][3]; //array of vectors of direct lattice A[vector number][coordinate(x,y,z) number]
        double B[3][3]; //array of vectors of reciprocal lattice B[vector number][coordinate(x,y,z) number]
        double ANORM[3]; //array of vector lengths (direct lattice) ANORM[vector number]
        double BNORM[3]; //array of vector lengths (reciprocal lattice) BNORM[vector number]
        double OMEGA; //cell volume (in direct space)
};

int latticeprint(latt);
int latticeclean(latt*);
int lattice_set_nan(latt*);
int latticecopy(latt, latt*);
int latticeprintfile(latt, ofstream*);

#endif
