#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include "vectors_fwd.h"

typedef struct
{
  int length;
  string_t string;
} string_appender;

void init_string_appender (string_appender *appender, string_t str);

void appender_strcat (string_appender *appender, const char *text);

#endif /* STRING_UTILS_H */
