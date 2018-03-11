#include "laspack_matrix.h"

int laspack_matrix_init (laspack_matrix *mat, const sparse_base_format *base)
{
  int row;
  int nnz = 0;

  Q_Constr (&mat->raw, "Matrix", base->N, False, Rowws, Normal, True);

  for (row = 0; row < base->N; row++)
    {
      int i;

      Q_SetLen (&mat->raw, row + 1, base->nnz_in_rows[row]);

      for (i = 0; i < base->nnz_in_rows[row]; i++)
        {
          Q_SetEntry (&mat->raw, row + 1, i, base->column_indecies[nnz + i] + 1, base->values[nnz + i]);
        }
      nnz += base->nnz_in_rows[row];
    }

  return 1;
}

void laspack_matrix_dump (laspack_matrix *mat, FILE *fout)
{
  int size = mat->raw.Dim;
  int i;
  for (i = 0; i < size; i++)
    {
      int j;
      for (j = 0; j < size; j++)
        {
          fprintf (fout, "%.3f ", Q_GetEl (&mat->raw, i + 1, j + 1));
        }
      fprintf (fout, "\n");
    }

  fflush (fout);
}

void laspack_matrix_destroy (laspack_matrix *mat)
{
  Q_Destr (&mat->raw);
}
