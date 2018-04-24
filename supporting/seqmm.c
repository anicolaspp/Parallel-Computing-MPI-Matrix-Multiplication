/*
 * Sequential program for matrix multiplication. Compile it using 
 * "gcc -o seqmm seqmm.c".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Define the following macro if you want to use the recursive block
   decomposition scheme for matrix multiplications; otherwise, it's
   going to be the dumb cache-unfriendly method. */
#define I_AM_SMART

typedef double datatype;

/* Cache size (in bytes): you need to modify this value to
   reflect the type of the machine on which you want to run this
   program. */
#define CACHE_SIZE 16384

/* read a matrix from a file */
void read_matrix(char* fname, datatype*** a, datatype** sa, int* m, int* n)
{
  FILE* finptr;
  int i, sz;

  finptr = fopen(fname, "r");
  if(!finptr) {
    perror("ERROR: can't open matrix file\n");
    exit(1); 
 }

  if(fread(m, sizeof(int), 1, finptr) != 1 ||
     fread(n, sizeof(int), 1, finptr) != 1) {
    perror("Error reading matrix file");
    exit(1);
  }
  sz = (*m)*(*n);

  *sa = (datatype*)malloc(sz*sizeof(datatype));
  if(fread(*sa, sizeof(datatype), sz, finptr) != sz) {
    perror("Error reading matrix file");
    exit(1);
  }

  *a = (datatype**)malloc((*m)*sizeof(datatype*));
  for(i=0; i<*m; i++) (*a)[i] = &(*sa)[i*(*n)];

  fclose(finptr);
}

/* write a matrix to a file */
void write_matrix(char* fname, datatype* sa, int m, int n)
{
  FILE* foutptr;
  int i;
  datatype* ptr;

  foutptr = fopen(fname, "w");
  if(!foutptr) {
    perror("ERROR: can't open matrix file\n");
    exit(1); 
 }

  if(fwrite(&m, sizeof(int), 1, foutptr) != 1 ||
     fwrite(&n, sizeof(int), 1, foutptr) != 1) {
    perror("Error reading matrix file");
    exit(1);
  }

  ptr = sa;
  for(i=0; i<m; i++) {
    if(fwrite(ptr, sizeof(datatype), n, foutptr) != n) {
      perror("Error reading matrix file");
      exit(1);
    }
    ptr += n;
  }

  fclose(foutptr);
}

/* dumb matrix multiplication; used for debugging purposes */
void dumb_matmul(datatype** a, datatype** b, datatype** c, int N) {
  int i, j, k;
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      for(k=0; k<N; k++)
	c[i][j] += a[i][k]*b[k][j];
}

/* matrix multiplication using recursive block decomposition (from Quinn's book) */
void recursive_matmul(int crow, int ccol, /* corner of C block */
	       int arow, int acol, /* corner of A block */
	       int brow, int bcol, /* corner of B block */
	       int l, int m, int n, /* block A is l*m, block B is m*n, block C is l*n */
	       int N, /* matrices are N*N */
	       datatype** a, datatype** b, datatype** c) {
  int i, j, k; 
  int lhalf[3], mhalf[3], nhalf[3]; /* quadrant sizes */
  datatype *aptr, *bptr, *cptr;

  if(m*n*sizeof(datatype) > CACHE_SIZE) { /* block B doesn't fit in cache */
    lhalf[0] = 0; lhalf[1] = l/2; lhalf[2] = l-l/2;
    mhalf[0] = 0; mhalf[1] = m/2; mhalf[2] = m-m/2;
    nhalf[0] = 0; nhalf[1] = n/2; nhalf[2] = n-n/2;
    for(i=0; i<2; i++)
      for(j=0; j<2; j++)
	for(k=0; k<2; k++)
	  recursive_matmul(crow+lhalf[i], ccol+nhalf[j],
			   arow+lhalf[i], acol+mhalf[k],
			   brow+mhalf[k], bcol+nhalf[j],
			   lhalf[i+1], mhalf[k+1], nhalf[j+1],
			   N, a, b, c);
  } else { /* block B fits in cache */
    for(i=0; i<l; i++) {
      for(j=0; j<n; j++) {
	cptr = &c[crow+i][ccol+j];
	aptr = &a[arow+i][acol];
	bptr = &b[brow][bcol+j];
	for(k=0; k<m; k++) {
	  *cptr += *(aptr++) * (*bptr);
	  bptr += N;
	}
      }
    }
  }
}

/* call this function to do matrix multiplication */
void matmul(datatype** a, datatype** b, datatype** c, int N) { 
#ifdef I_AM_SMART
  printf("use the recursive block decomposition scheme...\n");
  recursive_matmul(0, 0, 0, 0, 0, 0, N, N, N, N, a, b, c);
#else 
  printf("use the dumb method...\n");
  dumb_matmul(a, b, c, N);
#endif
}

int main (int argc, char * argv[]) 
{
  int n; /* dimension of the matrix */
  datatype *sa, *sb, *sc; /* storage for matrix A, B, and C */
  datatype **a, **b, **c; /* 2-d array to access matrix A, B, and C */
  int i, j;

  if(argc != 4) {
    printf("Usage: %s fileA fileB fileC\n", argv[0]);
    return 1;
  }

  /* read matrix A */
  printf("read matrix A from %s\n", argv[1]);
  read_matrix(argv[1], &a, &sa, &i, &j);
  if(i != j) { printf("ERROR: matrix A not square\n"); return 2; }
  n = i;

  /* read matrix B */
  printf("read matrix B from %s\n", argv[2]);
  read_matrix(argv[2], &b, &sb, &i, &j);
  if(i != j) { printf("ERROR: matrix B not square\n"); return 2; }
  if(n != i) { printf("ERROR: matrix A and B incompatible\n"); return 2; }

  /* initialize matrix C */
  sc = (datatype*)malloc(n*n*sizeof(datatype));
  memset(sc, 0, n*n*sizeof(datatype));
  c = (datatype**)malloc(n*sizeof(datatype*));
  for(i=0; i<n; i++) c[i] = &sc[i*n];

  /* do the multiplication */
  matmul(a, b, c, n);
  
  /* write matrix C */
  write_matrix(argv[3], sc, n, n);
  /*
  for (i=0; i<n; i++) {
    for(j=0; j<n; j++) printf("%4d ", c[i][j]);
    printf("\n");
  }
  */

  printf("matrix multiplication result stored in %s\n", argv[3]);
  return 0;
}
