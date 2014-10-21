//
//  pnicorr_debug.h
//  pni_correlation_service
//
//  Created by Benjamin Singer on 10/16/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#ifndef _pnicorr_debug_h
#define _pnicorr_debug_h

#ifndef DEBUG
#define NDEBUG
#endif

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

// ******* Debugging
#ifdef NDEBUG
// release version; make debug mode commands no-ops
#define NOP (void)0
#define LOG(...) NOP
#define LOGIF(x, ...) NOP
#define TIC NOP
#define TOC(desc) NOP
#else
// for debugging
#define LOG(...) printf(__VA_ARGS__)
#define LOGIF(x, ...)                                                          \
  if (x)                                                                       \
  printf(__VA_ARGS__)
// timing
static struct timeval time1, time2;
#define TIC gettimeofday(&time1, NULL)
#define TOC(desc)                                                              \
  do {                                                                         \
    gettimeofday(&time2, NULL);                                                \
    pnicorr_reporttime(time1, time2, desc);                                    \
  } while (0)
#endif

// for feedback in release mode as well
#define ERRIF(cond, ...)                                                       \
  do {                                                                         \
    if (cond) {                                                                \
      sprintf(errmsg, __VA_ARGS__);                                            \
      perror(errmsg);                                                          \
      pnicorr_debugbreak(errmsg, __FILE__, __LINE__);                          \
    }                                                                          \
  } while (0)
#define ALOG(...) printf(__VA_ARGS__)
#define ALOGIF(x, ...)                                                         \
  if (x)                                                                       \
  printf(__VA_ARGS__)

#define MAX_MSGLEN 1024

static char errmsg[MAX_MSGLEN];

void pnicorr_debugbreak(const char *s, const char *file, int line);
void pnicorr_reporttime(const struct timeval start, const struct timeval end,
                        const char *desc);
#endif
