#include "table_io_private.h"
#include "common/vectors.h"
#include "common/debug_utils.h"
#include <string.h>



void table_io_init (table_io *handle,
                    int rows,
                    int cols,
                    char *entries[],
                    table_output_format_t format)
{
  handle->format = format;
  handle->rows = rows;
  handle->cols = cols;
  handle->length = 0;
  handle->end_index = 0;

  table_io_fill_cols_max (handle, entries);

  table_io_allocate_string (handle);

  table_io_make_table (handle, entries);
}

int table_io_latex_additional_len (table_io *handle, int rows, int cols)
{
  int len = 0;
  FIX_UNUSED (handle);

  len += strlen ("begin{tabular}\n");
  len += strlen ("{}\n");
  len += cols; /*c*/
  len += cols + 1; /*|*/
  len += rows;
  len += strlen ("\\hline\n") * (rows + 1);
  len += strlen ("end{tabular}");

  return len;
}

void table_io_fill_cols_max (table_io *handle, char *entries[])
{
  int i;
  int j;
  int max = 0;

  handle->cols_max_entry_length = VECTOR_CREATE (int, handle->cols);

  for (j = 0; j < handle->cols; j++)
    {
      int max_len = 0;
      for (i = 0; i < handle->rows; i++)
        {
          int cur_len;
          const char *entry = entries[i * handle->cols + j];

          cur_len = strlen (entry);

          if (cur_len > max_len)
            max_len = cur_len;
        }
      handle->cols_max_entry_length[j] = max_len;
    }

  for (j = 0; j < handle->cols; j++)
    max = handle->cols_max_entry_length[j] > max ? handle->cols_max_entry_length[j] : max;

  handle->buf_text = VECTOR_CREATE (char, max + 2);
}

void table_io_allocate_string (table_io *handle)
{
  int j;
  int entries_in_row_length = 0;

  if (handle->format == latex_format)
    handle->length += table_io_latex_additional_len (handle, handle->rows, handle->cols);


  for (j = 0; j < handle->cols; j++)
    entries_in_row_length += handle->cols_max_entry_length[j] + 4;

  handle->length += entries_in_row_length * handle->rows;

  handle->length += handle->rows * (handle->cols + 2);

  handle->table_text = VECTOR_CREATE (char, handle->length);
}

void table_io_make_table (table_io *handle, char *entries[])
{
  int rows = handle->rows;
  int cols = handle->cols;
  int i;
  int j;

  table_io_put_prefix (handle);
  for (i = 0; i < rows; i++)
    {
      table_io_begin_newline (handle);
      table_io_put_beginline_separator (handle);
      for (j = 0; j < cols; j++)
        {
          table_io_put_entry (handle, entries[i * cols + j], j);
          if (j != cols - 1)
            table_io_put_separator (handle);
        }
      table_io_put_endline_separator (handle);
    }
  table_io_begin_newline (handle);
  table_io_put_postfix (handle);
}

void table_io_put_prefix (table_io *handle)
{
  int j;
  if (handle->format == latex_format)
    {
      table_io_strcat (handle, "begin{tabular}{|");
      for (j = 0; j < handle->cols; j++)
        table_io_strcat (handle, "c|");
      table_io_strcat (handle, "}");
    }
}

void table_io_put_postfix (table_io *handle)
{
  if (handle->format == latex_format)
    table_io_strcat (handle, "end{tabular}");
}

void table_io_begin_newline (table_io *handle)
{
  if (handle->format == latex_format)
    table_io_strcat (handle, "\n\\hline\n");
  else
    table_io_strcat (handle, "\n");
}

void table_io_put_beginline_separator (table_io *handle)
{
  if (handle->format == latex_format)
    table_io_strcat (handle, " ");
  if (handle->format == human_readable)
    table_io_strcat (handle, "|");
}

void table_io_put_separator (table_io *handle)
{
  if (handle->format == latex_format)
    table_io_strcat (handle, "&");
  if (handle->format == human_readable)
    table_io_strcat (handle, "|");
  if (handle->format == gnuplot_xyz)
    table_io_strcat (handle, " ");
}

void table_io_put_entry (table_io *handle, const char *entry, int col)
{
  int space_pref_len = 0;
  int space_post_len = 0;
  int entry_len = strlen (entry);
  int must_be_len = handle->cols_max_entry_length[col] + 2;
  int dif_len = must_be_len - entry_len;
  space_pref_len = dif_len / 2;
  space_post_len = dif_len - space_pref_len;

  VECTOR_SET (char, handle->buf_text, ' ',  space_pref_len);

  table_io_strcat (handle, handle->buf_text);

  table_io_strcat (handle, entry);

  VECTOR_SET (char, handle->buf_text, ' ',  space_post_len);

  table_io_strcat (handle, handle->buf_text);
}

void table_io_put_endline_separator (table_io *handle)
{
  if (handle->format == latex_format)
    table_io_strcat (handle, "\\\\");
  if (handle->format == human_readable)
    table_io_strcat (handle, "|");
}

void table_io_destroy (table_io *handle)
{
  VECTOR_DESTROY (handle->table_text);
  VECTOR_DESTROY (handle->buf_text);
  VECTOR_DESTROY (handle->cols_max_entry_length);
}

void table_io_strcat (table_io *handle, const char* entry)
{
  int i_e = 0;

  while (entry[i_e] != 0)
    {
      handle->table_text[handle->end_index] = entry[i_e];
      handle->end_index++;
      i_e++;
    }

  handle->table_text[handle->end_index] = 0;
}
