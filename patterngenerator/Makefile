COMPILER = gcc
DFLAGS =
CFLAGS = -c
OFLAGS = -lpthread -lm -lz -o
EXECNAME = patternconverter
EXECNAME2 = patterntestreader

all:   patternconverter patterntestreader

patternconverter:   main.o functions.o
	${COMPILER} ${DFLAGS} main.o functions.o ${OFLAGS} ${EXECNAME}

patterntestreader:   patterntestreader_main.o functions.o
	${COMPILER} ${DFLAGS} patterntestreader_main.o functions.o ${OFLAGS} ${EXECNAME2}

main.o:   main.c functions.h
	${COMPILER} ${DFLAGS} ${CFLAGS} main.c

patterntestreader_main.o:   patterntestreader_main.c functions.h
	${COMPILER} ${DFLAGS} ${CFLAGS} patterntestreader_main.c

functions.o:   functions.c functions.h
	${COMPILER} ${DFLAGS} ${CFLAGS} functions.c

