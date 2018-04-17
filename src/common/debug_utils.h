#ifndef DEBUG_UTILS_H
#define DEBUG_UTILS_H

#ifdef __GNUC__

  #include <signal.h>

  #ifdef DEBUG

    #define DEBUG_PAUSE(reason) (raise (SIGTRAP));

  #else

    #define DEBUG_PAUSE(reason);

  #endif

#else

  #ifdef __WIN32__

    #define DEBUG_PAUSE(reason) (__debugbreak ());

  #else

    #define DEBUG_PAUSE(reason);

  #endif /*__WIN32__ */

#endif /*__GNUC>__ */

#define DEBUG_ASSERT(assertion) if (!(assertion)) {DEBUG_PAUSE ("Assertion failed")}
#define ASSERT_RETURN(assertion, retval) if (!(assertion)) {DEBUG_PAUSE (""); return retval;}
#define ASSERT_RETURN_VOID(assertion) if (!(assertion)) {DEBUG_PAUSE ("Assertion failed"); return;}

#define FIX_UNUSED (void)
#endif /* DEBUG_UTILS_H */
