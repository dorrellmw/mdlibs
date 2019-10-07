#ifndef PSFPDB
#define PSFPDB

struct atom {
  char seg[11];
  char resID[9];
  char resType[9];
  char name[9];
  double charge;
};

struct psfpdb {
  int natom;
  struct atom * atoms;
};

struct psfpdb readStruct(const char * path);
void freeStruct(struct psfpdb p);
#endif
