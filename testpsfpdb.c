#include <stdio.h>
#include <limits.h>

#include "psfpdb.h"

int main(int argc, const char* argv[]) {
  struct psfpdb p = readStruct(argv[1]);

  if(p.natom > 0) {
    printf("%d atoms found.\n",p.natom);

    int atomnum = INT_MAX % p.natom;
    printf("Sample information for atom %d:\n", atomnum);
    struct atom * a = &(p.atoms[atomnum]);
    printf("  Segment name: %s\n",a->seg);
    printf("  Residue ID: %s\n", a->resID);
    printf("  Residue name: %s\n", a->resType);
    printf("  Atom name: %s\n",a->name);
    printf("  Charge: %lf\n", a->charge);
  } else {
    printf("No atoms found.\n");
  }

  freeStruct(p);
  return 0;
}
