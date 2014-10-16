
#include <stdio.h>
#include <signal.h>
#include "pnicorr_debugbreak.h"

void pnicorr_debugbreak(const char *s, const char *file, int line) {
  fprintf(stderr, "Assertion failure: %s, at %s:%d\n", s, file, line);
  fflush(stderr);
#if defined(WIN32)
  /*
   * We used to call DebugBreak() on Windows, but amazingly, it causes
   * the MSVS 2010 debugger not to be able to recover a call stack.
   */
  *((volatile int *)NULL) = 0;
  exit(3);
#elif defined(__APPLE__)
  /*
   * On Mac OS X, Breakpad ignores signals. Only real Mach exceptions are
   * trapped.
   */
  *((volatile int *)NULL) =
      0;          /* To continue from here in GDB: "return" then "continue". */
  raise(SIGABRT); /* In case above statement gets nixed by the optimizer. */
#else
  raise(SIGABRT); /* To continue from here in GDB: "signal 0". */
#endif
}
