// In this task, the object is to employ guassian elimination method to find the solution of linear
//system  equations of matrix form
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openblas/lapacke.h>
//#include <mkl_lapacke.h>

double *generate_matrix(int size)
{
  int i;
  double *matrix = (double *) malloc(sizeof(double) * size * size);

  srand(1);

  for (i = 0; i < size * size; i++) {
    matrix[i] = rand() % 100;
  }

  return matrix;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return abs(x - y) <= epsilon * abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

int check_result(double *bref, double *b, int size)
{
  int i;

  for(i = 0; i < size*size; i++) {
    if (!is_nearly_equal(bref[i], b[i]))
      return 0;
  }

  return 1;
}

int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{

//Replace next line to use your own DGESV implementation
//  LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);



//Start of my dgesv routine

//Declaration of variables and arrays

int i,j,c,r,rank;
double MatArray[n],IdentityArray[n],term;
double sum=0;

//dynamic allocation of memory for the matrix and identity matrix
       //Pointers

        double *ident,*matri;

         //allocate memory

        ident = (double *) malloc(sizeof(double) * n * n);
        matri = (double *) malloc(sizeof(double) * n * n);


//identity matrix such that identity=1

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i == j) {
                    ident[i*n+j] = 1;
                }
                else {
                        ident[i*n+j]=0;
                }
          }


        }


//Ensuring the diagonal of matrix equal to 1 and making the upper and lower corners to 0

    for (i = 0; i < n; i++){
        for (j = 0; j <= i; j++){
        if (i>j){

         term = a[i*n+j]/a[rank*n+j];
          for (r= 0; r < n; r++){
         MatArray[r] = a[rank*n+r]*term;
            a[i*n+r] = a[i*n+r] - MatArray[r];
         MatArray[r] = ident[rank*n+r]*term;
        ident[i*n+r] = ident[i*n+r] - IdentityArray[r];
                                    }
          }

          else if (i == j){
            double term;
            term = a[i*n+j];
            for(c = 0; c < n; c++){
                a[i*n+c] = a[i*n+c]/term;
                ident[i*n+c] = ident[i*n+c]/term;
                        }
                  }
      }
    }

    for(i = (n-1); i >= 0; i--){
        for (j = (n-1); j >= i; j--){
        if (i<j){


        term = a[i*n+j]/a[rank*n+j];
          for (r = (n-1); r >= 0; r--){
       MatArray[r] = a[rank*n+r]*term;
          a[i*n+r] = a[i*n+r] - MatArray[r];
        IdentityArray[r] = ident[rank*n+r]*term;
        ident[i*n+r] = ident[i*n+r] - IdentityArray[r];




                                  }
        }
      }
    }


 //finally the inverse of a x b
    for(i = 0; i < n; i++){
        for (j = 0; j < n; j++)
                              {
                    for(c = 0; c < n; c++){
            sum = sum +a[i*n+c] * b[c*n+j];
                    }
            matri[i*n+j] = sum;
            }
    }

    for (i = 0; i < n * n; i++)

{
      b[i] = matri[i];

}
  // End of my_dgesv  routine
}

void main(int argc, char *argv[])
{
  int size = atoi(argv[1]);

  double *a, *aref;
  double *b, *bref;

  a = generate_matrix(size);
  aref = generate_matrix(size);
  b = generate_matrix(size);
  bref = generate_matrix(size);

  // Using LAPACK dgesv OpenBLAS implementation to solve the system
  int n = size, nrhs = size, lda = size, ldb = size, info;
  int *ipiv = (int *) malloc(sizeof(int) * size);

  clock_t tStart = clock();
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
  printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  int *ipiv2 = (int *) malloc(sizeof(int) * size);

  tStart = clock();
  my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  if (check_result(bref, b, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");
}
