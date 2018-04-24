CC = mpicc

SRC = main.c 
OBJ = $(SRC:.c=.o)

all:	hw4

hw4:	$(OBJ)
	@echo LINK $(OBJ) INTO $@
	$(CC) $(OBJ) -o $@

clean:
	rm -f *.o *~ hw4
