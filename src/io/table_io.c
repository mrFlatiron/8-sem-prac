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
}

void table_io_allocate_string (table_io *handle)
{
  int j;
  int entries_in_row_length = 0;

  if (handle->format == latex_format)
    handle->length += table_io_latex_additional_len (handle, handle->rows, handle->cols);


  for (j = 0; j < handle->cols; j++)
    entries_in_row_length += handle->cols_max_entry_length[j] + 2;

  handle->length += entries_in_row_length * handle->rows;

  handle->length += handle->rows * (handle->cols + 1);

  handle->length += handle->rows;

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
      strcat (handle->table_text, "begin{tabular}{|");
      for (j = 0; j < handle->cols; j++)
        strcat (handle->table_text, "c|");
      strcat (handle->table_text, "}");
    }
}

void table_io_put_postfix (table_io *handle)
{
  if (handle->format == latex_format)
    strcat (handle->table_text, "end{tabular}");
}

void table_io_begin_newline (table_io *handle)
{
  if (handle->format == latex_format)
    strcat (handle->table_text, "\n\\hline\n");
  else
    strcat (handle->table_text, "\n");
}

void table_io_put_beginline_separator (table_io *handle)
{
  if (handle->format == latex_format)
    strcat (handle->table_text, " ");
  if (handle->format == human_readable)
    strcat (handle->table_text, "|");
}

void table_io_put_separator (table_io *handle)
{
  if (handle->format == latex_format)
    strcat (handle->table_text, "&");
  if (handle->format == human_readable)
    strcat (handle->table_text, "|");
  if (handle->format == gnuplot_xyz)
    strcat (handle->table_text, " ");
}

void table_io_put_entry (table_io *handle, const char *entry, int col)
{
  int i;
  int space_pref_len = 0;
  int space_post_len = 0;
  int entry_len = strlen (entry);
  int must_be_len = handle->cols_max_entry_length[col] + 2;
  int dif_len = must_be_len - entry_len;
  space_pref_len = dif_len / 2;
  space_post_len = dif_len - space_pref_len;

  for (i = 0; i < space_pref_len; i++)
    strcat (handle->table_text, " ");

  strcat (handle->table_text, entry);

  for (i = 0; i < space_post_len; i++)
    strcat (handle->table_text, " ");
}

void table_io_put_endline_separator (table_io *handle)
{
  if (handle->format == latex_format)
    strcat (handle->table_text, "\\\\");
  if (handle->format == human_readable)
    strcat (handle->table_text, "|");
}

void table_io_destroy (table_io *handle)
{
  VECTOR_DESTROY (handle->table_text);
}
