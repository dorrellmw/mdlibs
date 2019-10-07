#include <stdio.h>
#include <limits.h>

#include "pdb.h"


int main(int argc, const char* argv[]) {
  struct pdb p = readPDB(argv[1]);

  if(p.cell.valid) {
    printf("Unit cell information (%lfx%lfx%lf):\n",p.cell.a,p.cell.b,p.cell.c);
    printf("  alpha = %lf\n",p.cell.alpha);
    printf("  beta  = %lf\n",p.cell.beta);
    printf("  gamma = %lf\n",p.cell.gamma);
    printf("  Space group: %s\n",p.cell.sGroup);
    printf("  Z value = %d\n",p.cell.z);

    printf("\n");
  }
  if(p.natom > 0) {
    printf("%d atoms found.\n",p.natom);

    int atomnum = INT_MAX % p.natom;
    printf("Sample information for atom %d:\n", atomnum);
    struct atom * a = &(p.atoms[atomnum]);
    printf("  Standard serial number: %d\n",a->serial);
    printf("  Nonstandard serial string: %s\n", a->mserial);
    printf("  Atom name: %s\n",a->name);
    printf("  Alternate location indicator: %c\n", a->altLoc);
    printf("  Standard residue name: %s\n", a->resName);
    printf("  Nonstandard residue name: %s\n", a->mresName);
    printf("  Chain identifier: %c\n", a->chainID);
    printf("  Standard residue sequence number: %d\n", a->resSeq);
    printf("  Nonstandard residue sequence number: %d\n", a->mresSeq);
    printf("  Insertion code: %c\n", a->iCode);
    printf("  Atom location: (%lf,%lf,%lf)\n", a->x,a->y,a->z);
    printf("  Occupancy: %lf\n", a->occupancy);
    printf("  Temperature factor: %lf\n", a->tempFactor);
    printf("  Element: %s\n", a->element);
    printf("  Charge: %s\n", a->charge);
    printf("  Nonstandard segment name: %s\n",a->mseg);
  } else {
    printf("No atoms found.\n");
  }

  freePDB(p);
  return 0;
}
