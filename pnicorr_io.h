//
//  pnicorr_io.h
//  pnicorr
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#ifndef __h__pnicorr_io__
#define __h__pnicorr_io__

#define MAX_LINE 4194304L

typedef enum {
  pnicorr_iotype_1D,
  pnicorr_iotype_1Dgz,
  pnicorr_iotype_mat,
  pnicorr_iotype_numiotypes
} pnicorr_iotype_t;

float *pnicorr_load_1D(const char *filename, int *nts, int *ntrs);
void pnicorr_savematrix(const float *result, const int M, const int N,
                        const char *mode, const pnicorr_iotype_t iotype,
                        char *outfile);

#ifdef __linux__

#define pnicorr_fmemopen fmemopen
#define pnicorr_strchrnul strchrnul

#else

FILE *pnicorr_fmemopen(void *buf, size_t size, const char *mode);
char *pnicorr_strchrnul(const char *s, int c);

#endif // __linux__

#endif /* defined(__h__pnicorr_io__) */
