CFLAGS = -std=c99

all: testpsf testpdb

testpsf: testpsf.c psf.c

testpdb: testpdb.c pdb.c
