
#include <stdio.h>
#include <signal.h>
#include "pnicorr_debug.h"

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

// ****************** timediff ***************
static struct timeval timediff(const struct timeval start,
                               const struct timeval end) {
  struct timeval temp;
  if ((end.tv_usec - start.tv_usec) < 0) {
    temp.tv_sec = end.tv_sec - start.tv_sec - 1;
    temp.tv_usec = 1000000 + end.tv_usec - start.tv_usec;
  } else {
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_usec = end.tv_usec - start.tv_usec;
  }
  return temp;
}

void pnicorr_reporttime(const struct timeval start, const struct timeval end,
                        const char *desc) {
  const struct timeval diff = timediff(start, end);

  if (diff.tv_sec == 0) {
    if (diff.tv_usec < 1000) {
      ALOG("Time for %s : %d usec\n", desc, diff.tv_usec);
    } else {
      ALOG("Time for %s : %d msec\n", desc, diff.tv_usec / 1000);
    }
  } else {
    ALOG("Time for %s : %ld sec\n", desc, diff.tv_sec);
  }
}
