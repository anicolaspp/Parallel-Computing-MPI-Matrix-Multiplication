CFLAGS = -O2
MPICC = mpicc

.PHONY:	all clean

all: cannon genmat prtmat seqmm

cannon:	cannon.c
	$(MPICC) $(CFLAGS) -o $@ $< -lm

genmat: genmat.c
	$(CC) $(CFLAGS) -o $@ genmat.c

prtmat: prtmat.c
	$(CC) $(CFLAGS) -o $@ prtmat.c

seqmm: seqmm.c
	$(CC) $(CFLAGS) -o $@ seqmm.c

clean:
	$(RM) cannon genmat prtmat seqmm
	$(RM) *.o core *~
