//
//  pnicorr_io.c
//  pni_correlation_service
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <zlib.h>
#include <float.h>
#include <assert.h>

#include "pnicorr_debug.h"
#include "pnicorr_io.h"

// =============    helper functions =====================
#pragma mark - writematmatrix section

/*
 *  writematmatrix and supporting functions are adapted from "matfiles.c"
 *  by Malcolm MacLean. The original can be found on Matlab Central here:

 https://www.mathworks.com/matlabcentral/fileexchange/26731-portable-matfile-exporter--in-c-

*/

#define miINT8 1
#define miUINT16 4
#define miINT32 5
#define miUINT32 6
#define miDOUBLE 9
#define miMATRIX 14

#define mxCELL_CLASS 1
#define mxCHAR_CLASS 4
#define mxDOUBLE_CLASS 6

/*
 * write a double to a stream in ieee754 format regardless of host
 *  encoding.
 *  x - number to write
 *  fp - the stream
 *  bigendian - set to write big bytes first, elee write litle bytes
 *              first
 *  Returns: 0 or EOF on error
 *  Notes: different NaN types and negative zero not preserved.
 *         if the number is too big to represent it will become infinity
 *         if it is too small to represent it will become zero.
 */
static inline int fwriteieee754(double x, FILE *fp, int bigendian) {
  int shift;
  unsigned long sign, exp, hibits, hilong, lowlong;
  double fnorm, significand;
  int expbits = 11;
  int significandbits = 52;

  /* zero (can't handle signed zero) */
  if (x == 0) {
    hilong = 0;
    lowlong = 0;
    goto writedata;
  }
  /* infinity */
  if (x > DBL_MAX) {
    hilong = 1024 + ((1 << (expbits - 1)) - 1);
    hilong <<= (31 - expbits);
    lowlong = 0;
    goto writedata;
  }
  /* -infinity */
  if (x < -DBL_MAX) {
    hilong = 1024 + ((1 << (expbits - 1)) - 1);
    hilong <<= (31 - expbits);
    hilong |= (1 << 31);
    lowlong = 0;
    goto writedata;
  }
  /* NaN - dodgy because many compilers optimise out this test, but
   *there is no portable isnan() */
  if (x != x) {
    hilong = 1024 + ((1 << (expbits - 1)) - 1);
    hilong <<= (31 - expbits);
    lowlong = 1234;
    goto writedata;
  }

  /* get the sign */
  if (x < 0) {
    sign = 1;
    fnorm = -x;
  } else {
    sign = 0;
    fnorm = x;
  }

  /* get the normalized form of f and track the exponent */
  shift = 0;
  while (fnorm >= 2.0) {
    fnorm /= 2.0;
    shift++;
  }
  while (fnorm < 1.0) {
    fnorm *= 2.0;
    shift--;
  }

  /* check for denormalized numbers */
  if (shift < -1022) {
    while (shift < -1022) {
      fnorm /= 2.0;
      shift++;
    }
    shift = -1023;
  }
  /* out of range. Set to infinity */
  else if (shift > 1023) {
    hilong = 1024 + ((1 << (expbits - 1)) - 1);
    hilong <<= (31 - expbits);
    hilong |= (sign << 31);
    lowlong = 0;
    goto writedata;
  } else
    fnorm = fnorm - 1.0; /* take the significant bit off mantissa */

  /* calculate the integer form of the significand */
  /* hold it in a  double for now */

  significand = fnorm * ((1LL << significandbits) + 0.5f);

  /* get the biased exponent */
  exp = shift + ((1 << (expbits - 1)) - 1); /* shift + bias */

  /* put the data into two longs (for convenience) */
  hibits = (long)(significand / 4294967296);
  hilong = (sign << 31) | (exp << (31 - expbits)) | hibits;
  lowlong = (unsigned long)(significand - hibits * 4294967296);

writedata:
  /* write the bytes out to the stream */
  if (bigendian) {
    fputc((hilong >> 24) & 0xFF, fp);
    fputc((hilong >> 16) & 0xFF, fp);
    fputc((hilong >> 8) & 0xFF, fp);
    fputc(hilong & 0xFF, fp);

    fputc((lowlong >> 24) & 0xFF, fp);
    fputc((lowlong >> 16) & 0xFF, fp);
    fputc((lowlong >> 8) & 0xFF, fp);
    fputc(lowlong & 0xFF, fp);
  } else {
    fputc(lowlong & 0xFF, fp);
    fputc((lowlong >> 8) & 0xFF, fp);
    fputc((lowlong >> 16) & 0xFF, fp);
    fputc((lowlong >> 24) & 0xFF, fp);

    fputc(hilong & 0xFF, fp);
    fputc((hilong >> 8) & 0xFF, fp);
    fputc((hilong >> 16) & 0xFF, fp);
    fputc((hilong >> 24) & 0xFF, fp);
  }
  return ferror(fp);
}

/*
 * put a 32 bit little-endian integer to file
 **/
static int fput32le(FILE *fp, size_t x) {
  int i;

  for (i = 0; i < 4; i++) {
    fputc(x & 0xFF, fp);
    x >>= 8;
  }
  return 0;
}

/*
 * Quick and easy function to write a double array as a MATLAB matrix file.
 * Params: fname - name of file to write
 *         name - MATLAB name of variable
 *         data - the data (in C format)
 *         m - number of rows
 *         n - number of columns
 *         transpose - flag for transposing data (i.e. data is in MATLAB format)
 *  Returns: 0 on success, -1 on fail
 */
static int writematmatrix(const char *fname, const char *name,
                          const float *data, const int m, const int n,
                          const int transpose) {
  FILE *fp;
  int i, ii;
  char buff[128];
  size_t bufflen;
  size_t totalsize;
  int err;

  assert(m > 0);
  assert(n > 0);
  assert((m * n) / n == m);
  fp = fopen(fname, "wb");
  if (!fp)
    return -1;
  sprintf(buff, "MATLAB matrix file of %s, generated by Malcolm McLean", name);
  bufflen = strlen(buff);
  for (i = 0; i < 123; i++)
    fputc(i < (int)bufflen ? buff[i] : ' ', fp);
  fputc(0, fp);

  fputc(0, fp);
  fputc(1, fp);
  fputc('I', fp);
  fputc('M', fp);

  /* the main tag */
  totalsize = 8 + 8 + strlen(name) + (n * m * 8) + 8 * 4;
  if (strlen(name) % 8)
    totalsize += 8 - (strlen(name) % 8);
  fput32le(fp, miMATRIX);
  fput32le(fp, totalsize);

  /* array descriptor field */
  fput32le(fp, miUINT32);
  fput32le(fp, 8);
  fputc(mxDOUBLE_CLASS, fp);
  fputc(4, fp); /* array flags */
  fputc(0, fp);
  fputc(0, fp);
  fput32le(fp, 0);

  /* array dimensions */
  fput32le(fp, miINT32);
  fput32le(fp, 2 * 4);
  fput32le(fp, m);
  fput32le(fp, n);

  /* array name */
  fput32le(fp, miINT8);
  fput32le(fp, strlen(name));
  for (i = 0; name[i]; i++)
    fputc(name[i], fp);
  while (i % 8) {
    fputc(0, fp);
    i++;
  }

  /* the actual data */
  fput32le(fp, miDOUBLE);
  fput32le(fp, m * n * 8);
  if (transpose) {
    for (i = 0; i < m; i++)
      for (ii = 0; ii < n; ii++)
        fwriteieee754(data[i * n + ii], fp, 0);
  } else {
    for (i = 0; i < n; i++)
      for (ii = 0; ii < m; ii++)
        fwriteieee754(data[ii * n + i], fp, 0);
  }

  if (ferror(fp)) {
    fclose(fp);
    return -1;
  }

  err = fclose(fp);
  if (err)
    return -1;

  return 0;
}

#pragma mark - linux compatibility

/*

 ported from
 https://github.com/ingenuitas/python-tesseract/blob/master/fmemopen.c

 which is missing; now see
 https://github.com/NimbusKit/memorymapping/blob/master/src/fmemopen.c

 */
#ifndef __linux__

struct fmem {
  size_t pos;
  size_t size;
  char *buffer;
};
typedef struct fmem fmem_t;

static int readfn(void *handler, char *buf, int size) {
  int count = 0;
  fmem_t *mem = handler;
  int available = (int)(mem->size - mem->pos);

  if (size > available)
    size = available;
  for (count = 0; count < size; mem->pos++, count++)
    buf[count] = mem->buffer[mem->pos];

  return count;
}

static int writefn(void *handler, const char *buf, int size) {
  int count = 0;
  fmem_t *mem = handler;
  size_t available = mem->size - mem->pos;

  if (size > (int)available)
    size = (int)available;
  for (count = 0; count < size; mem->pos++, count++)
    mem->buffer[mem->pos] = buf[count];

  return count; // ? count : size;
}

static fpos_t seekfn(void *handler, fpos_t offset, int whence) {
  fpos_t pos;
  fmem_t *mem = handler;

  switch (whence) {
  case SEEK_SET:
    pos = offset;
    break;
  case SEEK_CUR:
    pos = mem->pos + offset;
    break;
  case SEEK_END:
    pos = mem->size + offset;
    break;
  default:
    return -1;
  }

  if (pos < 0 || (size_t)pos > mem->size)
    return -1;

  mem->pos = pos;
  return pos;
}

static int closefn(void *handler) {
  free(handler);
  return 0;
}

/* simple, but portable version of fmemopen for OS X / BSD */
static FILE *fmemopen(void *buf, size_t size, const char *mode) {
  fmem_t *mem = calloc(1, sizeof(fmem_t));
  mem->size = size, mem->buffer = buf;
  return funopen(mem, readfn, writefn, seekfn, closefn);
}

/* taviso / rarvmtools */
/*   strchrnul */
/*  for non-linux platforms */
static char *strchrnul(const char *s, int c) {
  char *ptr = strchr(s, c);
  if (!ptr) {
    ptr = strchr(s, '\0');
  }
  return ptr;
}

#endif // if not __linux__

#pragma mark - public exports

// ****************** get_uncompressed_size ***************
static size_t get_uncompressed_size(const int file, int *is_compressed) {
  size_t filesize;
  int dupfile = dup(file);
  ERRIF(dupfile == -1, "dup");

  gzFile gzfile = gzdopen(dupfile, "r");
  ERRIF(NULL == gzfile, "gzdopen");

  int is_direct = gzdirect(gzfile);

  FILE *fp = fdopen(dupfile, "r");
  ERRIF(NULL == fp, "fdopen");

  if (is_direct) {

    fseek(fp, 0, SEEK_END);
    filesize = ftell(fp);
    rewind(fp);

  } else {

    // gzseek(... SEEK_END) is not supported by zlib.
    // can get uncompressed size from final 4 bytes in file
    // *** warning : won't work if file > 4GB! But input data is rarely that
    // big. this fn only does input

    ERRIF(fseek(fp, -4, SEEK_END) != 0, "fseek");

    unsigned char bufSize[4];
    size_t len = fread(&bufSize[0], sizeof(unsigned char), 4, fp);

    ERRIF(len != 4, "gzipped file doesn't have uncompressed size (modulo 4GB) "
                    "encoded in final 4 bytes");

    rewind(fp);

    filesize = (size_t)((bufSize[3] << 24) | (bufSize[2] << 16) |
                        (bufSize[1] << 8) | bufSize[0]);
  }

  fclose(fp); // dup so original is still open
  *is_compressed = (is_direct == 0);
  return filesize;
}

// ****************** get_trs_from_line ***************
static int get_trs_from_line(const char *line) {
  int trs = 0;
  char *next;
  const char *remaining = line;
  const char *pos = strchr(remaining, '.');
  // <untested
  char c;
  while (pos > line) {
    c = *pos;
    if (c == ' ' || c == '\t')
      break;
    --pos;
  }
  remaining = pos;
  // /untested>
  float val = strtof(remaining, &next);
  while (val || (remaining != next)) {
    ++trs;
    remaining = next;
    val = strtof(remaining, &next);
  }
  return trs;
}

// warning: get_vals_from_line advances the input data pointer to the end of the
// line!
// ****************** get_vals_from_line ***************
static int get_vals_from_line(const char *line, const int start_tr,
                              float **data_in) {
  int vals = 0;
  int els = 0;
  float *dptr = NULL;
  if (data_in) {
    dptr = *data_in;
  }
  const char *remaining = line;
  char *next;
  float val = strtof(remaining, &next);
  while (val || (remaining != next)) {
    if (dptr) {
      if (vals >= start_tr) {
        *dptr = val;
        dptr++;
        ++els;
      }
    }
    remaining = next;
    val = strtof(remaining, &next);
    ++vals;
  }
  if (dptr) {
    *data_in = dptr;
    return els;
  }

  return vals;
}

// ****************** save_result ***************
static void write_1Dgz(const float *result, const int M, const int N,
                       const char *mode, const char *outfile) {
  char compress_mode[4];
  const int default_compression_level = 4;

  sprintf(compress_mode, "%s%d", mode, default_compression_level);

  gzFile *fp = gzopen(outfile, compress_mode);
  ERRIF(NULL == fp, "gzopen");
  int i, j, start;
  for (j = 0; j < M; ++j) {
    start = j * N;
    for (i = 0; i < N; ++i) {
      gzprintf(fp, "%4.4f ", result[start + i]);
    }
    gzprintf(fp, "\n");
  }
  gzclose(fp);
}

static void write_1D(const float *result, const int M, const int N,
                     const char *mode, const char *outfile) {
  FILE *fp = fopen(outfile, mode);
  ERRIF(NULL == fp, "fopen");
  int i, j, start;
  for (j = 0; j < M; ++j) {
    start = j * N;
    for (i = 0; i < N; ++i) {
      fprintf(fp, "%4.4f ", result[start + i]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

float *pnicorr_load_1D(const char *filename, int *nts, int *ntrs) {
  // ************** read file into memory *****************
  int file = open(filename, O_RDONLY);
  ERRIF(-1 == file, "open");

  int is_compressed;
  size_t filesize = get_uncompressed_size(file, &is_compressed);

  LOG("file is %zu bytes in size\n", filesize);

  char *filebuf = (char *)malloc(filesize);
  ERRIF(NULL == filebuf, "malloc");

  gzFile gzfile = gzdopen(file, "r");
  ERRIF(NULL == gzfile, "gzdopen");
  size_t bytes_read = gzread(gzfile, filebuf, (int)filesize);
  ERRIF(bytes_read != filesize, "gzread");

  gzclose(gzfile); // also closes 'file'
  LOG("File read\n");

  // get the number of newlines in file
  // will be the number of num_timeseries + comment lines at top of file
  int nlines = 0;
  char *fb = filebuf;
  char *pos;
  size_t remaining = filesize;
  while ((pos = memchr(fb, '\n', remaining)) != NULL) {
    ++nlines;
    remaining -= (pos - fb + 1);
    fb = pos + 1;
    LOGIF((0 == (nlines % 10000)), "%d\n", nlines);
  }
  LOG("Counted %d lines\n", nlines);

  // ********************* parse data into time-series ***********
  char line[MAX_LINE];
  FILE *sp = fmemopen(filebuf, filesize, "r");
  ERRIF(NULL == sp, "fmemopen");

  int comment_lines = 0;
  char *next;
  float val;
  while (fgets(line, MAX_LINE, sp) != NULL) {
    val = strtof(line, &next);
    if (val || (line != next))
      break;
    ++comment_lines;
  }
  LOG("Counted %d comment lines\n", comment_lines);
  int num_timeseries = nlines - comment_lines;
  LOG("Subtracting comment lines, that's %d num_timeseries\n", num_timeseries);

  // Determine number of TRs in first non-comment line
  int vals = get_vals_from_line(line, 0, NULL);
  int trs = get_trs_from_line(line);
  int start_tr = vals - trs;
  LOG("Starting at TR %d, skipping prior columns\n", start_tr);
  LOG("Counted %d floating point valued TRs in time series\n", trs);

  ALOG("matrix is %d num_timeseries x %d trs\nLoading ...\n", num_timeseries,
       trs);

  // Knowing num_timeseries and trs, we can allocate memory for the float data
  int n_elements = num_timeseries * trs;
  size_t datasize = n_elements * sizeof(float);
  LOG("Allocating %zu bytes to hold the data\n", datasize);
  assert(datasize > 0);
  float *data = (float *)malloc(datasize);
  ERRIF(NULL == data, "malloc");

  LOG("Data allocated.\n");

  LOG("Reading %d floats...\n", n_elements);
  int i;
  int el;
  const int update_count =
      (num_timeseries > 100) ? (num_timeseries / 100) : num_timeseries;
  int count = 0;
  float *dptr = data;
  for (i = 0; i < num_timeseries; ++i) {
    el = get_vals_from_line(line, start_tr, &dptr);
    // dptr now points to the next row-- on final iteration
    // this means it is garbage!
    if (el != trs) {
      LOG("%s\n", line);
      LOG("%d\n", el);
    }
    ERRIF(el != trs, "Data for node %d incomplete\n", i);
    count += el;
    ALOGIF(0 == (i % update_count), "%d%%\n",
           (int)((float)i / (float)num_timeseries * 100.0f));
    fgets(line, MAX_LINE, sp);
  }
  LOG("first 10:\n");
  assert(data != NULL);
  for (i = 0; i < 10; ++i) {
    LOG("%4f\n", data[i]);
  }
  LOG("last 10:\n");
  for (i = count - 10; i < count; ++i) {
    LOG("%4f\n", data[i]);
  }
  ERRIF(count != n_elements, "something wrong: read %d elements, but it should "
                             "have been %d x %d or %d elements",
        count, trs, num_timeseries, n_elements);

  fclose(sp);
  free(filebuf);

  *nts = num_timeseries;
  *ntrs = trs;
  return data;
}

const char *pnicorr_genoutfile(const char *filename, const int num_timeseries,
                               const char *ext, const char *comp) {
  char *outfile = calloc(MAX_FILENAME, 1);
  char *loc = strchrnul(filename, '.');
  if (*loc != '\0')
    loc = strchr(filename, '.'); // want first one, if it's there at all
  if (loc[1] == '/')
    loc = strchr(loc, '.'); // but if the form ./file.1D, want second one
  char tmp[MAX_FILENAME];
  strncpy(&tmp[0], filename, loc - filename);
  sprintf(outfile, "%s_%dx%d_correlations.%s%s", tmp, num_timeseries,
          num_timeseries, ext, comp);
  return outfile;
}

// ****************** save_result ***************
void pnicorr_savematrix(const float *result, const int M, const int N,
                        const char *mode, const pnicorr_iotype_t iotype,
                        const char *outfile) {
  switch (iotype) {
  case pnicorr_iotype_1Dgz:
    write_1Dgz(result, M, N, mode, outfile);
    break;
  case pnicorr_iotype_1D:
    write_1D(result, M, N, mode, outfile);
    break;
  case pnicorr_iotype_mat:
    ERRIF(mode[0] == 'a', "cannot append to mat files!\n");
    writematmatrix(outfile, "pnicorr", result, M, N, 0);
    break;
  default:
    ERRIF(1, "unknown  iotype %d\n", iotype);
    break;
  }
}
