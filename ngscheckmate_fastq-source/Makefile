CC ?= gcc
DFLAGS = 
CFLAGS = -c
OFLAGS = -lpthread -lm -lz -o


EXECNAME = ngscheckmate_fastq

ngscheckmate:   stringhash2.o ngscheckmate_functions.o ngscheckmate_main.o
	${CC} ${DFLAGS} stringhash2.o ngscheckmate_functions.o ngscheckmate_main.o ${OFLAGS} ${EXECNAME}

stringhash2.o:   stringhash2.c stringhash2.h
	${CC} ${DFLAGS} ${CFLAGS} stringhash2.c

ngscheckmate_functions.o:   ngscheckmate_functions.c ngscheckmate.h patternconvert.h
	${CC} ${DFLAGS} ${CFLAGS} ngscheckmate_functions.c patternconvert.h

ngscheckmate_main.o:     ngscheckmate_main.c ngscheckmate.h
	${CC} ${DFLAGS} ${CFLAGS} ngscheckmate_main.c


