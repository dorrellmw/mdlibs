#ifndef PDB
#define PDB

struct pdbatom {
  int serial;
  char name[5];
  char altLoc;
  char resName[4];
  char chainID;
  int resSeq;
  char iCode;
  double x;
  double y;
  double z;
  double occupancy;
  double tempFactor;
  char element[3];
  char charge[3];
  //BONUS: nonstandard PDB field modifications
  char mserial[7];
  char mresName[5];
  int mresSeq;
  char mseg[11];
};

struct cryst {
  int valid;
  double a;
  double b;
  double c;
  double alpha;
  double beta;
  double gamma;
  char sGroup[12];
  int z;
};

struct pdb {
  struct cryst cell;
  int natom;
  struct pdbatom * atoms;
};

struct pdb readPDB(const char * path);
void freePDB(struct pdb p);
#endif
