# Define required macros here
SHELL = /bin/sh

EXE = audio_omp
OBJS =  ${EXE}.o
CC = gcc
CFLAGS = -w -O3 -fopenmp -march=armv8-a -mtune=cortex-a57
LDFLAGS =  -fopenmp
INCLUDE =
LIBS = 

${EXE}:${OBJS}
	${CC} ${LDFLAGS} ${INCLUDES} -o $@ ${OBJS} ${LIBS}

clean:
	-rm -f *.o core *.core ${EXE}

.c.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<
