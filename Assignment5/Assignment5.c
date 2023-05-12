#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_lapacke.h>

void read_input(long *N_x, long *N_y, double *L_x, double *L_y, double *t_f, double *t_D, int *s_x0, int *s_x1, int *s_y0, int *s_y1);

void read_coefficients(long N_x, long N_y, double *u);

long indx(long j, long p, long J);

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

long N_x; //number of x grid points.//
long N_y; //number of y grid points.//
double L_x; //Length of x domain.//
double L_y; //Length of y domain.//
double t_f; //final time for time evolution mode.//
double t_D; //diagnostic timestep.//
int s_x0; //switch for x = 0 boundary condition.//
int s_x1; //switch for x = Lx boundary condition.//
int s_y0; //switch for y = 0 boundary condition.//
int s_y1; //switch for y = Ly boundary condition.//


    /* Read in from file; */
read_input(&N_x, &N_y, &L_x, &L_y, &t_f, &t_D, &s_x0, &s_x1, &s_y0, &s_y1);

/* Grid spacing */
double dx;
dx = L_x/(N_x-1);

double dy;
dy = L_y/(N_y-1);

double invdx2 = 1.0/(dx*dx);

double invdy2 = 1.0/(dy*dy);

double *u = (double*)malloc(sizeof(double)*(N_x*N_y));

if (u == NULL) {
        printf("Failed to allocated %lu bytes memory", (N_x*N_y));
        exit(1);
}

read_coefficients(N_x, N_y, u);

//I have no idea what im doing//

   band_mat bmat;
   long nbands_low = N_x;
   long nbands_up = N_x;
   init_band_mat(&bmat, nbands_low, nbands_up, (N_x*N_y));
   double *x = malloc(sizeof(double)*(N_x*N_y));
   double *b = malloc(sizeof(double)*(N_x*N_y));



   

/* Loop over the equation number and set the matrix values
equal to the coefficients of the grid values note
boundaries treated with special cases. */

long j,p;

int deez;

    for(j=0;j<N_x;j++) {
        for(p=0;p<N_y;p++) {
            if(j==0 && p==0) {
                deez = s_x0*s_y0;
            } else if(j==0 && p==N_y-1) {
                deez = s_x0*s_y1;
            } else if(j==N_x-1 && p==0) {
                deez = s_x1*s_y0;
            } else if(j==N_x-1 && p==N_y-1) {
                deez = s_x1*s_y1;
            } else if(j==0) {
                deez = s_x0;
            } else if(j==N_x-1) {
                deez = s_x1;
            } else if(p==0) {
                deez = s_y0;
            } else if(p==N_y-1) {
                deez = s_y1;
            } else {
                deez = 1;
            }
            setv(&bmat,indx(j,p,N_x),indx(j,p,N_x),1+2*t_D*invdx2+2*t_D*invdy2);
            if(j==0) {
                setv(&bmat, indx(j,p,N_x), indx(j+1,p,N_x) , -deez*(t_D*invdx2));
                setv(&bmat, indx(j,p,N_x), indx(j,p,N_x), deez*(1+t_D*invdx2+2*t_D*invdy2));
            } 
            else if(j<N_x-1) {
                setv(&bmat,indx(j,p,N_x),indx(j-1,p,N_x),-t_D*invdx2);
                setv(&bmat,indx(j,p,N_x),indx(j+1,p,N_x),-t_D*invdx2);
            }
            else {
                setv(&bmat, indx(j,p,N_x), indx(j-1,p,N_x), -deez*(t_D*invdx2));
                setv(&bmat, indx(j,p,N_x), indx(j,p,N_x), deez*(1+t_D*invdx2+2*t_D*invdy2));
            }
            if(p==0) {
                setv(&bmat, indx(j,p,N_x), indx(j,p+1,N_x), -deez*(t_D*invdy2));
                setv(&bmat, indx(j,p,N_x), indx(j,p,N_x), deez*(1+2*t_D*invdx2+t_D*invdy2));
            } 
            else if(p<N_y-1) {
                setv(&bmat,indx(j,p,N_x),indx(j,p-1,N_x),-t_D*invdy2);
                setv(&bmat,indx(j,p,N_x),indx(j,p+1,N_x),-t_D*invdy2);
            }
            else {
                setv(&bmat, indx(j,p,N_x), indx(j,p-1,N_x), -deez*(t_D*invdy2));
                setv(&bmat, indx(j,p,N_x), indx(j,p,N_x), deez*(1+2*t_D*invdx2+t_D*invdy2));
            }
            if(indx(j,p,N_x) == 0 || indx(j,p,N_x) == N_x-1 || indx(j,p,N_x) == N_x*N_y-N_x || indx(j,p,N_x) == N_x*N_y-1) {
                setv(&bmat, indx(j,p,N_x), indx(j,p,N_x), deez*(1+t_D*invdx2+t_D*invdy2));
            }
            if(getv(&bmat, indx(j,p,N_x), indx(j,p,N_x))==0) {
                setv(&bmat, indx(j,p,N_x), indx(j,p,N_x), 1);
            }

            if(j==0) {
                u[indx(j,p,N_x)] = s_x0*u[indx(j,p,N_x)];
            }
            if(p==0) {
                u[indx(j,p,N_x)] = s_y0*u[indx(j,p,N_x)];
            }
            if(j==N_x-1) {
                u[indx(j,p,N_x)] = s_x1*u[indx(j,p,N_x)];
            }
            if(p==N_y-1) {
                u[indx(j,p,N_x)] = s_y1*u[indx(j,p,N_x)];
            }

        }
    }

            /* source term */
        for(long k=0; k<N_x*N_y; k++){
        b[k] = u[k] + (t_D*(u[k]*u[k]))/(1+(u[k]*u[k]));
        x[k] = u[k];
        }




FILE *output_file;
    if(!(output_file=fopen("output.txt","w"))) {
       printf("Error creating an output\n");
       exit(1);
    }
double *tmp;
long t;

    for(t=0;t*t_D<t_f;t++) {

        for(long p=0;p<N_y;p++) {
            for(long j=0;j<N_x;j++) {
                fprintf(output_file, "%lf %lf %lf %lf\n", t*t_D, j*dx, p*dy, x[indx(j,p,N_x)]);
            }
        }
        solve_Ax_eq_b(&bmat, x, b);
        tmp = x;
        x = b;
        b = tmp;
        for(long j=0;j<N_x*N_y;j++) {
            b[j] = b[j] + t_D*b[j]*b[j]/(1+b[j]*b[j]);
        }
    }

    fclose(output_file);


}

void read_input(long *N_x, long *N_y, double *L_x, double *L_y, double *t_f, double *t_D, int *s_x0, int *s_x1, int *s_y0, int *s_y1) {
   FILE *input_file;
   if(!(input_file=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(10!=fscanf(input_file,"%ld %ld %lf %lf %lf %lf %i %i %i %i", N_x, N_y, L_x, L_y, t_f, t_D, s_x0, s_x1, s_y0, s_y1)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(input_file);
}

void read_coefficients(long N_x, long N_y, double *u) {
        FILE *coefficients_file;
    if(!(coefficients_file=fopen("coefficients.txt", "r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    for (long i = 0; i < N_x*N_y; i++) {
        if(1!=fscanf(coefficients_file, "%lf", &u[i])) {
            printf("Error reading coefficients from file\n");
            exit(1);
        }
    }
    fclose(coefficients_file);
}

long indx(long j, long p, long J) {
  return p*J + j;
}
