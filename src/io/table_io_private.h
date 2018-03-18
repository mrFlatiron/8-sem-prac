#ifndef TABLE_IO_PRIVATE_H
#define TABLE_IO_PRIVATE_H

#include "table_io.h"

void table_io_fill_cols_max (table_io *handle, char *entries[]);
void table_io_allocate_string (table_io *handle);
void table_io_make_table (table_io *handle, char *entries[]);

int table_io_latex_additional_len (table_io *handle, int rows, int cols);

void table_io_put_prefix (table_io *handle);
void table_io_put_postfix (table_io *handle);
void table_io_begin_newline (table_io *handle);
void table_io_put_beginline_separator (table_io *handle);
void table_io_put_separator (table_io *handle);
void table_io_put_entry (table_io *handle, const char *entry, int col);
void table_io_put_endline_separator (table_io *handle);

#endif /* TABLE_IO_PRIVATE_H */
