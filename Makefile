CFLAGS = -std=c99

all: testpsf testpdb testpsfpdb

testpsf: testpsf.c psf.c

testpdb: testpdb.c pdb.c

testpsfpdb: testpsfpdb.c psfpdb.c psf.c pdb.c

.PHONY: clean
clean:
	-rm -f testpsfpdb testpsf testpdb
