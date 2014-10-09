#include "gsl_linalg.h"
#include "gsl_check_range.h"
#include "gsl_permute_vector.h"
#include "gsl_blas.h"

#include <math.h>

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno)
{
	fprintf(stderr, "Error %s (errno %d) in file %s, line %d", reason, gsl_errno, file, line);
	exit(1);
}

#define REAL double

#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)

#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir

#include "block/init_source.c"

#include "permutation/init.c"
#include "permutation/permutation.c"
#include "permutation/permute_source.c"

#include "matrix/init_source.c"
#include "matrix/matrix_source.c"
#include "matrix/copy_source.c"
#include "matrix/swap_source.c"

#include "vector/init_source.c"
#include "vector/vector_source.c"
#include "vector/copy_source.c"

int gsl_check_range;

int
gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int *signum)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("LU decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != A->size1)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      size_t i, j, k;

      *signum = 1;
      gsl_permutation_init (p);

      for (j = 0; j < N - 1; j++)
        {
          /* Find maximum in the j-th column */

          REAL ajj, max = fabs (gsl_matrix_get (A, j, j));
          size_t i_pivot = j;

          for (i = j + 1; i < N; i++)
            {
              REAL aij = fabs (gsl_matrix_get (A, i, j));

              if (aij > max)
                {
                  max = aij;
                  i_pivot = i;
                }
            }

          if (i_pivot != j)
            {
              gsl_matrix_swap_rows (A, j, i_pivot);
              gsl_permutation_swap (p, j, i_pivot);
              *signum = -(*signum);
            }

          ajj = gsl_matrix_get (A, j, j);

          if (ajj != 0.0)
            {
              for (i = j + 1; i < N; i++)
                {
                  REAL aij = gsl_matrix_get (A, i, j) / ajj;
                  gsl_matrix_set (A, i, j, aij);

                  for (k = j + 1; k < N; k++)
                    {
                      REAL aik = gsl_matrix_get (A, i, k);
                      REAL ajk = gsl_matrix_get (A, j, k);
                      gsl_matrix_set (A, i, k, aik - aij * ajk);
                    }
                }
            }
        }
      
      return GSL_SUCCESS;
    }
}

int
gsl_linalg_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LU->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for x */

      gsl_linalg_LU_svx (LU, p, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_LU_svx (const gsl_matrix * LU, const gsl_permutation * p, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution/rhs size", GSL_EBADLEN);
    }
  else
    {
      /* Apply permutation to RHS */

      gsl_permute_vector (p, x);

      /* Solve for c using forward-substitution, L c = P b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasUnit, LU, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LU, x);

      return GSL_SUCCESS;
    }
}
