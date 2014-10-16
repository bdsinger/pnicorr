//
//  pnicorr_reporttime.h
//  pnicorr
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#ifndef __h__pnicorr_reporttime__
#define __h__pnicorr_reporttime__

#include <sys/time.h>
static struct timeval time1, time2;
#define TIC gettimeofday(&time1, NULL)
#define TOC(desc)                                                              \
  do {                                                                         \
    gettimeofday(&time2, NULL);                                                \
    pnicorr_reporttime(time1, time2, desc);                                    \
  } while (0)

void pnicorr_reporttime(const struct timeval start, const struct timeval end,
                        const char *desc);

#endif /* defined(_h__pnicorr_reporttime__) */
