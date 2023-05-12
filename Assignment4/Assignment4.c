#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_lapacke.h>

void read_input(double *L, long *N, double *D, double *v, double *k_pos, double *k_neg);

void read_coefficients(long N, double *S, double *sigma);

long indexing(long i, long N);

void write_output(long N, double dx, double *A, double *B);

// Band matrix//

struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat */
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }  
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) { 
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}

int main() {

  /******************************/
  /* Declarations of parameters */
  /******************************/                 
  double L;
  long N;
  double D;
  double v;
  double k_pos;
  double k_neg;

  double *A, *B;

  /* Read in from file; */
  read_input(&L, &N, &D, &v, &k_pos, &k_neg);

  A = (double*)malloc(sizeof(double)*N);
  B = (double*)malloc(sizeof(double)*N);

 /* Grid spacing */
  double dx;
  if (N > 2){
    dx = L/N;
  }
  else {
    printf("Need a larger N. \n");
    return 1;
  }

  double invdx2 = 1.0/(dx*dx);


   double *S = (double*)malloc(sizeof(double)*N);
   double *sigma = (double*)malloc(sizeof(double)*N);

   if (S == NULL || sigma == NULL) {
        printf("Failed to allocated %lu bytes memory", 2*N);
        exit(1);
   }

   read_coefficients(N, S, sigma);

   //Set up matrix and account for boundary conditions//

   double alpha = D*invdx2;
   double beta_A = -(v/dx + k_pos + (2*D)*invdx2);
   double beta_B = -(v/dx + k_neg + (2*D)*invdx2);
   double gamma = (D*invdx2 + v/dx);

   band_mat bmat;
   long nbands_low = 4;
   long nbands_up = 4;
   init_band_mat(&bmat, nbands_low, nbands_up, 2*N);
   double *x = malloc(sizeof(double)*(2*N));
   double *b = malloc(sizeof(double)*(2*N));
   


   long i;
   for (i = 0; i<2*N; i++) {
    if(i>2*N-3) {setv(&bmat, indexing(i, 2*N), indexing(i-(2*N-2), 2*N), alpha);};
    if(i>1)       {setv(&bmat,indexing(i, 2*N), indexing(i-2, 2*N), gamma);};
    if(i<2*N-2) {setv(&bmat, indexing(i, 2*N), indexing(i+2, 2*N), alpha);};
    if(i<2)     {setv(&bmat, indexing(i, 2*N), indexing(i+(2*N-2), 2*N), gamma);};
   }

   for (i = 0; i < 2*N; i++) {
   if (i % 2 == 0) {setv(&bmat, indexing(i, 2*N), indexing(i, 2*N), beta_A);}
   else {setv(&bmat, indexing(i, 2*N), indexing(i, 2*N), beta_B - sigma[(i-1)/2]);};
   }
   
   for (i = 1; i < 2*N; i++) {
   if ((i-1) % 2 == 0) {setv(&bmat, indexing(i, 2*N), indexing(i-1, 2*N), k_pos);}
   else {setv(&bmat, indexing(i, 2*N), indexing(i-1, 2*N), 0);};
   }

   for (i = 0; i < 2*N-1; i++) {
   if (i % 2 == 0) {setv(&bmat, indexing(i, 2*N), indexing(i+1, 2*N), k_neg);}
   else {setv(&bmat, indexing(i, 2*N), indexing(i+1, 2*N), 0);};
   }

   for (i=0; i < 2*N; i++) {
     if (i % 2 == 0) {
       b[indexing(i, 2*N)] = -S[i/2];
     }
     else {
       b[indexing(i, 2*N)] = 0;
     }
   }

   //Solve multilinear equation//
   
   solve_Ax_eq_b(&bmat, x, b);
   
   for (i=0; i < N/2; i++) {
    A[i] = x[4*i];
    B[i] = x[(4*i)+2];
    }

   for (i = N/2; i < N; i++) {
   A[i] = x[(2*N-1)-4*(i-(N/2))];
   B[i] = x[((2*N-1)-4*(i-(N/2)))-2];
   }
   
   write_output(N, dx, A, B);

   finalise_band_mat(&bmat);


   free(x);
   free(b);
   free(S);
   free(sigma);
   free(A);
   free(B);
}

void read_input(double *L, long *N, double *D, double *v, double *k_pos, double *k_neg) {
   FILE *input_file;
   if(!(input_file=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(6!=fscanf(input_file,"%lf %ld %lf %lf %lf %lf", L, N, D, v, k_pos, k_neg)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(input_file);
}

 void read_coefficients(long N, double *S, double *sigma) {
    FILE *coefficients_file;
    if(!(coefficients_file=fopen("coefficients.txt", "r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    for (long i = 0; i < N; i++) {
        if(2!=fscanf(coefficients_file, "%lf %lf", &S[i], &sigma[i])) {
            printf("Error reading coefficients from file\n");
            exit(1);
        }
    }
    fclose(coefficients_file);
}

   long indexing(long i, long N) {
      if (i < N/2) {
        return 2*i;
        }
        else {
        return 2*(N-i)-1;
        }
      }

void write_output(long N, double dx, double *A, double *B) {
    FILE *output_file;
    if(!(output_file=fopen("output.txt","w"))) {
       printf("Error creating an output\n");
       exit(1);
    }
    for (long i = 0; i < N; i++) {
       fprintf(output_file, "%lf %lf %lf\n", dx*i, A[i], B[i]);
    }
   fclose(output_file);
  }

