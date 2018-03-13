#ifndef TABLE_IO_H
#define TABLE_IO_H

#include "kernel/kernel_typedefs.h"
#include "common/vectors_fwd.h"

/*
human readable:

| A\B | hor header | another header |
| vh  |     1      |      2.5       |
| vh2 |    234     |     123.45     |


latex:

\begin{center}
\begin{tabular}{|c|c|c|c|c|}
\hline
N\textbackslash M     &      20      &      40      &      80      &     160      \\
\hline
           20         & 3,922676e-03 & 6,846325e-03 & 7,855319e-03 & 8,028590e-03 \\
\hline
           80         & 2,284488e-03 & 6,666070e-04 & 1,385733e-03 & 1,552270e-03 \\
\hline
           320        & 3,298528e-03 & 5,208735e-04 & 1,538522e-04 & 3,205057e-04 \\
\hline
          1280        & 3,528191e-03 & 7,919148e-04 & 1,299347e-04 & 3,777771e-05 \\
\hline
\end{tabular}
\end{center}

 */

typedef struct
{
  string_t table_text;

  int length;

  int rows;
  int cols;

  vector_int_t cols_max_entry_length;

  table_output_format_t format;
} table_io;


void table_io_init (table_io *handle, int rows, int cols,
                    const char *entries[],
                    table_output_format_t format);


void table_io_destroy (table_io *handle);

#endif /* TABLE_IO_H */
