//
//  pnicorr_io.h
//  pnicorr
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//

#ifndef __h__pnicorr_io__
#define __h__pnicorr_io__

float *pnicorr_load_1D(const char *filename, int *nts, int *ntrs);
void pnicorr_savematrix(const float *result, const int M, const int N,
                        const char *mode, const int save_text,
                        const int compression_level, char *outfile);

#endif /* defined(__h__pnicorr_io__) */
