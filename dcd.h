#ifndef DCD_H_
#define DCD_H_

#include <stdio.h>
#include <stdint.h>

struct DCD;

struct DCD *openDCD(char *);
struct DCD *openWritableDCD(char *);
void closeDCD(struct DCD *);
uint32_t getNFrames(struct DCD *);
void goToFrame(struct DCD *,uint32_t);
void nextFrame(struct DCD *);
uint32_t getFrame(struct DCD *);
void getUnitCell(struct DCD *, double *);
void getCoords(struct DCD *, float *, float *, float *);
void writeCoords(struct DCD *, float *, float *, float *);

#endif
