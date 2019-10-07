#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "pdb.h"

/**
 * Read crystallographic information from a PDB
 *
 * Reads unit cell details from a PDB file and stores the contents in the
 * provided pdb struct.
 *
 * @param[in] pdb The PDB from which the unit cell data will be read.
 * @param[in,out] p The struct into which the unit cell data will be stored.
 * @return 0 if cell data is read successfully, or -1 if an error occurs.
 */
int readCryst(FILE * pdb, struct pdb * p) {
  rewind(pdb); // Start at the beginning of the file.
  char buffer[128]; // For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  while (1) {
    ptr = fgets(buffer,128,pdb);
    if(ferror(pdb) || feof(pdb) || ptr==NULL) {
      break; // Error reading file for cryst record.
    }
    if(buffer[0]=='\n') {
      continue; // Skip empty lines.
    }
    if(strncmp(buffer,"CRYST1",6)) {
      continue; // Skip non-CRYST1 lines.
    }

    int offset=0; // Keeps track of field alignment
    char substr[15];


    int width=6; // Width of first field ("CRYST1")
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the CRYST1 prefix
    if(strncmp(buffer,"CRYST1",6))
      continue; // Prefix doesn't match (this shouldn't ever happen)
    offset+=width; // advance offset to next field

    width=9; // Width of second field (a)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the 'a' substring
    p->cell.a=atof(substr); // write 'a' value to cryst struct
    offset+=width; // advance offset to next field

    width=9; // Width of third field (b)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the 'b' substring
    p->cell.b=atof(substr); // write 'b' value to cryst struct
    offset+=width; // advance offset to next field

    width=9; // Width of fourth field (c)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the 'c' substring
    p->cell.c=atof(substr); // write 'c' value to cryst struct
    offset+=width; // advance offset to next field

    width=7; // Width of fifth field (alpha)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the alpha substring
    p->cell.alpha=atof(substr); // write alpha value to cryst struct
    offset+=width; // advance offset to next field

    width=7; // Width of sixth field (beta)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the beta substring
    p->cell.beta=atof(substr); // write beta value to cryst struct
    offset+=width; // advance offset to next field

    width=7; // Width of seventh field (gamma)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the gamma substring
    p->cell.gamma=atof(substr); // write gamma value to cryst struct
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=11; // Width of eighth field (space group)
    // copy string directly into cryst struct
    memcpy(p->cell.sGroup,&buffer[offset],width);
    offset+=width; // advance offset to next field

    width=4; // Width of ninth field (Z value)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the Z value substring
    p->cell.z=atoi(substr); // write resSeq value to cryst struct
    offset+=width; // advance offset to next field

    p->cell.valid=true; // Mark cryst struct as valid
    break; // Only process the first CRYST1 line, ignore any others
  } // This loop ends upon interruption by a break statement.
  if(p->cell.valid)
    return 0;
  return -1;
}

/**
 * Read atom information from a PDB
 *
 * Reads atom details from a PDB file and stores the contents in the provided
 * pdb struct. Also populates the natom field of the pdb struct.
 *
 * @param[in] pdb The PDB from which the atom data will be read.
 * @param[in,out] p The struct into which the atom data will be stored.
 * @return The number of atoms read
 */
int readAtoms(FILE * pdb, struct pdb * p) {
  rewind(pdb); // Start at the beginning of the file.
  char buffer[128]; // For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Allocate some space for atom array
  int size = 100;
  p->atoms = malloc(size * sizeof(struct pdbatom));

  int n = 0; // Track the number of atoms read so far

  while (1) {
    ptr = fgets(buffer,128,pdb);
    if(ferror(pdb) || feof(pdb) || ptr==NULL) {
      break; // Error reading file for next atom record.
    }
    if(buffer[0]=='\n') {
      continue; // Skip empty lines.
    }
    if(strncmp(buffer,"ATOM  ",6)) {
      continue; // Skip non-ATOM lines.
    }

    // Expand array if necessary
    if(n==size) {
      size*=2;
      p->atoms = realloc(p->atoms,size * sizeof(struct pdbatom));
    }

    // Clear pdbatom struct
    p->atoms[n].serial = -1;
    memset(p->atoms[n].name,'\0',5);
    p->atoms[n].altLoc = '\0';
    memset(p->atoms[n].resName,'\0',4);
    p->atoms[n].chainID = '\0';
    p->atoms[n].resSeq = 0;
    p->atoms[n].iCode = '\0';
    p->atoms[n].x = 0;
    p->atoms[n].y = 0;
    p->atoms[n].z = 0;
    p->atoms[n].occupancy = 0;
    p->atoms[n].tempFactor = 0;
    memset(p->atoms[n].element,'\0',3);
    memset(p->atoms[n].charge,'\0',3);
    memset(p->atoms[n].mserial,'\0',7);
    memset(p->atoms[n].mresName,'\0',5);
    p->atoms[n].mresSeq = 0;
    memset(p->atoms[n].mseg,'\0',11);

    int offset=0; // Keeps track of field alignment
    char substr[15];

    int width=6; // Width of first field ("ATOM  ")
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the ATOM prefix
    if(strncmp(buffer,"ATOM  ",6))
      continue; // Prefix doesn't match (this shouldn't ever happen)
    offset+=width; // advance offset to next field

    width=5; // Width of second field (serial)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the serial substring
    p->atoms[n].serial=atoi(substr); // write serial value to pdbatom struct
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=4; // Width of third field (name)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].name,&buffer[offset],width);
    offset+=width; // advance offset to next field

    width=1; // Width of fourth field (altLoc)
    // copy char directly into pdbatom struct
    p->atoms[n].altLoc = buffer[offset];
    offset+=width; // advance offset to next field

    width=3; // Width of fifth field (resName)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].resName,&buffer[offset],width);
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=1; // Width of sixth field (chainID)
    // copy char directly into pdbatom struct
    p->atoms[n].chainID = buffer[offset];
    offset+=width; // advance offset to next field

    width=4; // Width of seventh field (resSeq)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the resSeq substring
    p->atoms[n].resSeq=atoi(substr); // write resSeq value to pdbatom struct
    offset+=width; // advance offset to next field

    width=1; // Width of eighth field (iCode)
    // copy char directly into pdbatom struct
    p->atoms[n].iCode = buffer[offset];
    offset+=width; // advance offset to next field

    offset+=3; // Skip three spaces

    width=8; // Width of ninth field (x position)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the x position substring
    p->atoms[n].x=atof(substr); // write x value to pdbatom struct
    offset+=width; // advance offset to next field

    width=8; // Width of tenth field (y position)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the y position substring
    p->atoms[n].y=atof(substr); // write y value to pdbatom struct
    offset+=width; // advance offset to next field

    width=8; // Width of eleventh field (z position)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the z position substring
    p->atoms[n].z=atof(substr); // write z value to pdbatom struct
    offset+=width; // advance offset to next field

    width=6; // Width of twelfth field (occupancy)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the occupancy substring
    // write occupancy value to pdbatom struct
    p->atoms[n].occupancy=atof(substr);
    offset+=width; // advance offset to next field

    width=6; // Width of thirteenth field (tempFactor)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the tempFactor substring
    // write tempFactor value to pdbatom struct
    p->atoms[n].tempFactor=atof(substr);
    offset+=width; // advance offset to next field

    offset+=10; // Skip ten spaces

    width=2; // Width of fourteenth field (element)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].element,&buffer[offset],width);
    offset+=width; // advance offset to next field

    width=2; // Width of fifth field (charge)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].charge,&buffer[offset],width);
    offset+=width; // advance offset to next field

    //BONUS: Non-standard PDB modifications

    width=6; // Width of modified second field (serial)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].mserial,&buffer[6],width);

    width=4; // Width of modified fifth field (resName)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].mresName,&buffer[17],width);

    width=5; // Width of modified seventh field (resSeq)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[22],width); // isolate the resSeq substring
    p->atoms[n].mresSeq=atoi(substr); // write resSeq value to pdbatom struct

    width=10; // Width of extra field in blank span (segment name)
    // copy string directly into pdbatom struct
    memcpy(p->atoms[n].mseg,&buffer[66],width);

    n++; // Advance atom count
  } // This loop ends upon interruption by a break statement.
  p->atoms = realloc(p->atoms,n * sizeof(struct pdbatom)); // Trim array size
  p->natom = n; // Store atom count in struct
  return n;
}

/**
 * Read data from PDB
 *
 * Parses a PDB file and creates a pdb struct containing the successfully parsed
 * information, including crystallographic data, the array of atom records, and
 * the total atom count.
 *
 * In the event of a read error, the relevant information will be ignored, but
 * other data will continue to be read and stored in the pdb struct if possible.
 *
 * @param[in] path The filesystem path to the PDB to be parsed.
 * @return A pdb struct containing all of the parsed information.
 */
struct pdb readPDB(const char * path) {

  struct pdb p = { .cell = { .valid = false,
                             .a = 0,
                             .b = 0,
                             .c = 0,
                             .alpha = 0,
                             .beta = 0,
                             .gamma = 0,
                             .sGroup = "\0\0\0\0\0\0\0\0\0\0\0\0",
                             .z = -1 },
                   .natom = -1,
                   .atoms = NULL};

  FILE * pdb = fopen(path,"r");
  if(!pdb) // Error encountered while opening file.
    return p;

  // A proper PDB file should not have lines more than 80 chars wide.
  char buffer[128]; //For storing one line of the file at a time.

  char * ptr; // For detecting read failures

  ptr = fgets(buffer,128,pdb); //get first line
  if(ferror(pdb) || feof(pdb) || ptr==NULL) {
    fclose(pdb);
    return p; // Error reading first line of file.
  }

  // Read crystal data
  readCryst(pdb, &p);

  // Read atom data
  readAtoms(pdb, &p);

  fclose(pdb);
  return p;
}

/**
 * Frees the memory allocated for the pdb struct
 *
 * Checks to make sure memory has not already been freed, and if not, frees the
 * array used to store atom information. Once freed, it sets the pointer to the
 * allocation to NULL to prevent double freeing.
 *
 * @param[in] p the pdb struct to be freed.
 */
void freePDB(struct pdb p) {
  void sfree(void ** ptr) {
    if(*ptr) {
      free(*ptr);
      *ptr = NULL;
    }
  }
  sfree((void **) &p.atoms);
}
