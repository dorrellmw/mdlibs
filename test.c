#include <stdio.h>
#include <limits.h>

#include "psf.h"


int main(int argc, const char* argv[]) {
  struct psf p = readPSF(argv[1]);

  if(p.sig.valid)
    printf("PSF opened successfully.\n");
  else {
    printf("Error encountered while reading PSF header.\n");
    return -1;
  }

  printf("PSF extended? %s\n",p.sig.ext?"yes":"no");
  printf("CMAP CHEQ? %s\n",p.sig.cmapcheq?"yes":"no");
  printf("XPLOR? %s\n",p.sig.xplor?"yes":"no");
  printf("Scattering lengths present? %s\n",p.sig.slb?"yes":"no");

  printf("\n");

  if(p.ntitle != -1) {
    printf("Titles (%d):\n", p.ntitle);
    for (int i=0; i<p.ntitle; i++)
      printf("%s\n",p.titles[i]);

    printf("\n");
  }

  printf("Atom information %s\n",p.natom==-1?"not found":"found");
  if(p.natom!=-1) {

    printf("Number of atoms: %d\n",p.natom);
  
    int atomnum = INT_MAX % p.natom;
    printf("Sample information for atom %d:\n", atomnum);
    struct atom * a = &(p.atoms[atomnum]);
    printf("  Segment name: %s\n",a->seg);
    printf("  Residue number: %s\n",a->resid);
    printf("  Residue name: %s\n",a->res);
    printf("  Atom name: %s\n",a->name);
    printf("  Atom type: %s\n",a->type);
    printf("  Atom charge: %lf\n",a->charge);
    printf("  Atom mass: %lf\n",a->mass);
    switch(a->imove) {
      case 0:
        printf("  IMOVE parameter: free\n");
        break;
      case -1:
        printf("  IMOVE parameter: lone pair\n");
        break;
      case 1:
        printf("  IMOVE parameter: fixed\n");
        break;
      default:
        printf("  IMOVE parameter: unknown (%d)\n",a->imove);
        break;
    }
    if(p.sig.cmapcheq) {
      printf("  CHEQ electronegativity: %lf\n",a->ech);
      printf("  CHEQ hardness: %lf\n",a->eha);
    }
    if(p.sig.slb) {
      printf("  Atom scattering length: %lf\n",a->b);
    }

    printf("\n");

    printf("Bond information %s\n",p.nbond==-1?"not found":"found");
    if(p.nbond!=-1) {
      printf("Number of bonds: %d\n",p.nbond);
      printf("Sample bonds for atom %d:\n", atomnum);
      for(int i=0; i<p.nbond; i++) {
        if(p.bonds[i].a == atomnum || p.bonds[i].b == atomnum)
          printf("  %d -- %d\n", p.bonds[i].a, p.bonds[i].b);
      }

      printf("\n");
    }

    printf("Angle information %s\n",p.ntheta==-1?"not found":"found");
    if(p.ntheta!=-1) {
      printf("Number of angles: %d\n",p.ntheta);
      printf("Sample angles for atom %d:\n", atomnum);
      for(int i=0; i<p.ntheta; i++) {
        if(p.angles[i].a == atomnum || p.angles[i].b == atomnum || p.angles[i].c == atomnum)
          printf("  %d -- %d -- %d\n", p.angles[i].a,p.angles[i].b,p.angles[i].c);
      }

      printf("\n");
    }

    printf("Dihedral information %s\n",p.nphi==-1?"not found":"found");
    if(p.nphi!=-1) {
      printf("Number of dihedrals: %d\n",p.nphi);
      printf("Sample dihedrals for atom %d:\n", atomnum);
      for(int i=0; i<p.nphi; i++) {
        if(p.dihedrals[i].a == atomnum || p.dihedrals[i].b == atomnum ||
            p.dihedrals[i].c == atomnum || p.dihedrals[i].d == atomnum)
          printf("  %d -- %d -- %d -- %d\n", p.dihedrals[i].a,p.dihedrals[i].b,p.dihedrals[i].c,p.dihedrals[i].d);
      }

      printf("\n");
    }

    printf("Improper dihedral information %s\n",p.nimphi==-1?"not found":"found");
    if(p.nimphi!=-1) {
      printf("Number of impropers: %d\n",p.nimphi);
      printf("Sample impropers for atom %d:\n", atomnum);
      for(int i=0; i<p.nimphi; i++) {
        if(p.impropers[i].a == atomnum || p.impropers[i].b == atomnum ||
            p.impropers[i].c == atomnum || p.impropers[i].d == atomnum)
          printf("  %d -- %d -- %d -- %d\n", p.impropers[i].a,p.impropers[i].b,p.impropers[i].c,p.impropers[i].d);
      }
    }
  }

  freePSF(p);
  return 0;
}
