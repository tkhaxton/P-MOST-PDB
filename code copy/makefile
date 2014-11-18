
COMPILER=gcc
CFLAGS=-lm -O3
#CFLAGS=-lm -g

all : roundtrip create_molecular_files

roundtrip : roundtrip.c mystdlib.c mystdlib.h mymath.c mymath.h protein_functions.c protein_functions.h
	${COMPILER} ${CFLAGS} -o roundtrip roundtrip.c protein_functions.c mymath.c mystdlib.c

create_molecular_files : create_molecular_files.c mystdlib.c mystdlib.h mymath.c mymath.h protein_functions.c protein_functions.h
	${COMPILER} ${CFLAGS} -o create_molecular_files create_molecular_files.c protein_functions.c mymath.c mystdlib.c

clean:
	rm create_molecular_files
	rm roundtrip

