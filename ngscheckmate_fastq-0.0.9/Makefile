COMPILER = gcc
DFLAGS = 
CFLAGS = -c
OFLAGS = -lpthread -lm -lz -o


EXECNAME = ngscheckmate_fastq

ngscheckmate:   stringhash2.o ngscheckmate_functions.o ngscheckmate_main.o
	${COMPILER} ${DFLAGS} stringhash2.o ngscheckmate_functions.o ngscheckmate_main.o ${OFLAGS} ${EXECNAME}

stringhash2.o:   stringhash2.c stringhash2.h
	${COMPILER} ${DFLAGS} ${CFLAGS} stringhash2.c

ngscheckmate_functions.o:   ngscheckmate_functions.c ngscheckmate.h patternconvert.h
	${COMPILER} ${DFLAGS} ${CFLAGS} ngscheckmate_functions.c patternconvert.h

ngscheckmate_main.o:     ngscheckmate_main.c ngscheckmate.h
	${COMPILER} ${DFLAGS} ${CFLAGS} ngscheckmate_main.c


