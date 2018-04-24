/*
 * Print out the matrix stored in a file.  Simply compile it using
 * "gcc -o prtmat prtmat.c".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double datatype;

int main (int argc, char * argv[]) 
{
  int i, j;
  int m, n;
  FILE *finptr;
  datatype *a;
  int blocksize;
  
  /* the user can optionally choose to print out only the first few
     rows and columns of the matrix */
  if(argc != 2 && argc != 3) {
    printf("Usage: %s filename [blocksize]\n", argv[0]);
    return 1;
  }

  finptr = fopen (argv[1], "r");
  if(!finptr) {
    perror("ERROR: can't open matrix file\n");
    return 2;
  }

  if(argc == 3) {
    blocksize = atoi(argv[2]);
    if(blocksize < 0) {
      fprintf(stderr, "ERROR: invalid blocksize argument\n");
      return 2;
    }
  } else blocksize = 1000000; /* basically infinite */

  if(fread(&m, sizeof(int), 1, finptr) != 1 ||
     fread(&n, sizeof(int), 1, finptr) != 1) {
    perror("Error reading matrix file");
    return 3;
  }
  printf("Matrix[%d][%d]:\n", m, n);

  a = (datatype*)malloc(n*sizeof(datatype));
  for (i = 0; i < m && i < blocksize; i++) {
    if(fread(a, sizeof(datatype), n, finptr) != n) {
      perror("Error reading matrix file");
      return 3;
    }
    /*printf("[row=%3d] ", i);*/
    for(j=0; j<n && j<blocksize; j++) printf("%.2lf ", (double)a[j]);
    printf("\n");
  }

  free(a);
  fclose (finptr);
  return 0;
}
