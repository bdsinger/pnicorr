//
//  pnicorr_time.c
//  pni_correlation_service
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#include "pnicorr.h"
#include "pnicorr_reporttime.h"

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
      LOG("Time for %s : %d usec\n", desc, diff.tv_usec);
    } else {
      LOG("Time for %s : %d msec\n", desc, diff.tv_usec / 1000);
    }
  } else {
    LOG("Time for %s : %ld sec\n", desc, diff.tv_sec);
  }
}
