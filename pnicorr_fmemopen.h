//
//  pnicorr_fmemopen.h
//  pnicorr
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#ifndef _pnicorr_fmemopen_h
#define _pnicorr_fmemopen_h
#ifdef __cplusplus
#include <cstdio>
extern "C" {
#else
#include <stdio.h>
#endif

#ifdef __linux__

#define pnicorr_fmemopen fmemopen
#define pnicorr_strchrnul strchrnul

#else

FILE *pnicorr_fmemopen(void *buf, size_t size, const char *mode);
char *pnicorr_strchrnul(const char *s, int c);

#endif
#ifdef __cplusplus
}
#endif
#endif
