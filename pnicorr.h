//
//  pnicorr.h
//  pni_correlation_service
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#ifndef pni_correlation_service_pnicorr_h
#define pni_correlation_service_pnicorr_h

#include <stdio.h>
#include <stdlib.h>

#include "pnicorr_debugbreak.h"
#include "pnicorr_reporttime.h"

// ******* Memory limits
#define BYTES_PER_GB (1000000000LL)
#define BYTES_PER_MB (1000000LL)
#define MAX_DIGITS 32
#define MAX_FILENAME 256
#define MAX_MSGLEN 1024
#define MAX_LINE 4194304L
#define DEFAULT_GZIP_LEVEL 4

// ******* Debugging
#ifdef NDEBUG
// release version; make debug mode commands no-ops
#define LOG(...) (__ASSERT_VOID_CAST(0))
#define LOGIF(x, ...) (__ASSERT_VOID_CAST(0))
#define ERRIF(cond, ...) (__ASSERT_VOID_CAST(0))
#define TIC (__ASSERT_VOID_CAST(0))
#define TOC(desc) (__ASSERT_VOID_CAST(0))
#else
// for debugging
#define LOG(...) printf(__VA_ARGS__)
#define LOGIF(x, ...)                                                          \
  if (x)                                                                       \
  printf(__VA_ARGS__)
#define ERRIF(cond, ...)                                                       \
  do {                                                                         \
    if (cond) {                                                                \
      sprintf(errmsg, __VA_ARGS__);                                            \
      perror(errmsg);                                                          \
      pnicorr_debugbreak(errmsg, __FILE__, __LINE__);                          \
    }                                                                          \
  } while (0)
#endif

// for feedback in release mode as well
#define ALOG(...) printf(__VA_ARGS__)
#define ALOGIF(x, ...)                                                         \
  if (x)                                                                       \
  printf(__VA_ARGS__)

static char errmsg[MAX_MSGLEN];

// ******* Runtime options , enums
typedef enum {
  pnicorr_rflags_nonorm,
  pnicorr_rflags_savetext,
  pnicorr_rflags_compressoutput,
  pnicorr_rflags_memory,
  pnicorr_rflags_numflags
} pnicorr_rflags_t;

static inline const int *pnicorr_runtimeflags(void) {
  int *flags = malloc(pnicorr_rflags_numflags * sizeof(int));

  flags[pnicorr_rflags_nonorm] = 1 << 0;
  flags[pnicorr_rflags_savetext] = 1 << 1;
  flags[pnicorr_rflags_compressoutput] = 1 << 2;
  flags[pnicorr_rflags_memory] = 1 >> 4;

  return flags;
}

static inline const char **pnicorr_flagtext(void) {
  const char **text = malloc(pnicorr_rflags_numflags * sizeof(const char *));

  text[pnicorr_rflags_nonorm] = "-nonorm";
  text[pnicorr_rflags_savetext] = "-savetext";
  text[pnicorr_rflags_compressoutput] = "-gzout";
  text[pnicorr_rflags_memory] = "-mem";

  return text;
}

static inline const char **pnicorr_flagdescription(void) {
  const char **desc = malloc(pnicorr_rflags_numflags * sizeof(const char *));

  desc[pnicorr_rflags_nonorm] = "do not normalize rows";
  desc[pnicorr_rflags_savetext] = "save text rather than binary files";
  desc[pnicorr_rflags_compressoutput] =
      "compress the output\n\t\t-gzout1 .. -gzout9 sets "
      "level\n\t\t-gzout alone is same as -gzout4";
  desc[pnicorr_rflags_memory] =
      "memory (MB)\n\t\t-mem1 .. -mem999999 asks for 1M through "
      "nearly 1T.\n\t\tSmaller means more file activity; "
      "computing is done \n\t\tin stages. Default is the same as "
      "-mem4000  (4G)\n";

  return desc;
}

#endif
