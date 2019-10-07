#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdb.h"
#include "psf.h"
#include "psfpdb.h"

/**
 * Extracts psfpdb struct information from a PSF file
 *
 * Reads a psf file and extracts the basic atom information, ignoring other
 * kinds of information that might also be present in the PSF. If an error
 * occurs, nothing is written.
 *
 * @param[in] path The filesystem path to the PSF to be parsed.
 * @param[out] result The psfpdb struct into which the atom data is written
 */
static void processPSF(const char * path, struct psfpdb * result) {
  struct psf p = readPSF(path);
  if(!p.sig.valid)
    return;
  if(p.natom == -1)
    return;
  result->natom = p.natom;
  result->atoms = malloc(result->natom * sizeof(struct atom));
  for(int i=0; i<p.natom; i++) {
    memset(result->atoms[i].seg,'\0',11);
    memset(result->atoms[i].resID,'\0',9);
    memset(result->atoms[i].resType,'\0',9);
    memset(result->atoms[i].name,'\0',9);
    result->atoms[i].charge=0;

    strcpy(result->atoms[i].seg,p.atoms[i].seg);
    strcpy(result->atoms[i].resID,p.atoms[i].resid);
    strcpy(result->atoms[i].resType,p.atoms[i].res);
    strcpy(result->atoms[i].name,p.atoms[i].name);
    result->atoms[i].charge = p.atoms[i].charge;
  }
  freePSF(p);
  return;
}

/**
 * Extracts psfpdb struct information from a PDB file
 *
 * Reads a pdb file and extracts the basic atom information, ignoring other
 * kinds of information that might also be present in the PDB. If an error
 * occurs, nothing is written.
 *
 * @param[in] path The filesystem path to the PDB to be parsed.
 * @param[out] result The psfpdb struct into which the atom data is written
 */
static void processPDB(const char * path, struct psfpdb * result) {
  struct pdb p = readPDB(path);
  if(p.natom == -1)
    return;
  result->natom = p.natom;
  result->atoms = malloc(result->natom * sizeof(struct atom));
  for(int i=0; i<p.natom; i++) {
    memset(result->atoms[i].seg,'\0',11);
    memset(result->atoms[i].resID,'\0',9);
    memset(result->atoms[i].resType,'\0',9);
    memset(result->atoms[i].name,'\0',9);
    result->atoms[i].charge=0;

    strcpy(result->atoms[i].seg,p.atoms[i].mseg);
    sprintf(result->atoms[i].resID,"%d",p.atoms[i].mresSeq);
    strcpy(result->atoms[i].resType,p.atoms[i].mresName);
    strcpy(result->atoms[i].name,p.atoms[i].name);
    int abscharge = p.atoms[i].charge[0]-'0';
    if((0 <= abscharge && abscharge <= 9) && 
        (p.atoms[i].charge[1]=='-' || p.atoms[i].charge[1]=='+')) {
      result->atoms[i].charge = abscharge;
      if(p.atoms[i].charge[1]=='-')
        result->atoms[i].charge*=-1;
    }
  }
  freePDB(p);
  return;
}

/**
 * Reads atomic structure information from a PSF or PDB file.
 * 
 * Infers the file type from the path given, and uses the appropriate PSF or
 * PDB header to parse the file into a simple psfpdb struct containing basic
 * atom information.
 * 
 * Returns an empty psfpdb struct with natom set to -1 if an error occurs.
 *
 * @param[in] path The filesystem path to the PSF or PDB to be parsed.
 * @return A psfpdb struct containing basic atom information.
 */
struct psfpdb readStruct(const char * path) {
  char buf[5]="\0\0\0\0\0";
  // copy file extension from argument into buffer
  strcpy(buf,&path[strlen(path)-4]);
  struct psfpdb result = { .natom = -1, .atoms = NULL};
  if(!strcmp(buf,".psf")) {
    processPSF(path,&result);
  }
  else if(!strcmp(buf,".pdb")) {
    processPDB(path,&result);
  }
  return result;
}

/**
 * Frees the memory allocated for the structure
 *
 * Checks to make sure memory has not already been freed, and if not, frees the
 * array used to store atom information. Once freed, it sets the pointer to the
 * allocation to NULL to prevent double freeing.
 *
 * @param[in] p the psfpdb struct to be freed.
 */
void freeStruct(struct psfpdb p) {
  if(p.atoms) {
    free(p.atoms);
    p.atoms = NULL;
  }
}
