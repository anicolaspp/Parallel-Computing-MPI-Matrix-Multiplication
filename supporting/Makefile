CFLAGS = -O2
MPICC = mpicc

.PHONY:	all clean

all: genmat prtmat seqmm

genmat: genmat.c
	$(CC) $(CFLAGS) -o $@ genmat.c

prtmat: prtmat.c
	$(CC) $(CFLAGS) -o $@ prtmat.c

seqmm: seqmm.c
	$(CC) $(CFLAGS) -o $@ seqmm.c

clean:
	$(RM) genmat prtmat seqmm
	$(RM) *.o core *~
