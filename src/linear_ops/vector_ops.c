#include "vector_ops.h"
#include "common/vectors.h"
#include  <math.h>

double dot_product (vector_double_t a, vector_double_t b, int size)
{
  int i;
  double sum = 0;

  for (i = 0; i < size; i++)
    sum += a[i] * b[i];

  return sum;
}

void linear_combination_1 (const vector_double_t a,
                           double coef, const vector_double_t b,
                           vector_double_t out,
                           int size)
{
  int i;

  VECTOR_SET (double, out, 0, size);

  for (i = 0; i < size; i++)
    out[i] = a[i] + coef * b[i];
}



void linear_combination_w_override_1 (vector_double_t a, double coef, const vector_double_t b, int size)
{
  int i;

  for (i = 0; i < size; i++)
    a[i] += coef * b[i];
}

double l2_norm (const vector_double_t a, int size)
{
  int i;
  double sum = 0;

  for (i = 0; i < size; i++)
    sum += a[i] * a[i];

  return sqrt (sum);
}

double c_norm (const vector_double_t a, int size)
{
  int i;
  double max = 0;

  for (i = 0; i < size; i++)
    max = (max < fabs (a[i])) ? fabs (a[i]) : max;

  return max;
}

double c_norm_w_index (const vector_double_t a, int size, int *index)
{
  int i;
  double max = 0;

  for (i = 0; i < size; i++)
    {
      if (max < fabs (a[i]))
        {
          if (index)
            *index = i;

          max = fabs (a[i]);
        }
    }

  return max;
}
