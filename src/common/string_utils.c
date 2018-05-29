#include "string_utils.h"
#include "common/debug_utils.h"
#include <string.h>

void init_string_appender (string_appender *appender, string_t str)
{
  ASSERT_RETURN_VOID (str && appender);

  appender->length = strlen (str);
  appender->string = str;
}

void appender_strcat (string_appender *appender, const char *text)
{
  int i = 0;
  int j = appender->length;

  ASSERT_RETURN_VOID (text);

  while (text[i] != 0)
    {
      appender->string[j] = text[i];
      i++;
      j++;
    }
  appender->string[j] = 0;
  appender->length = j;
}
