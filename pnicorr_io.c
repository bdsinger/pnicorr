//
//  pnicorr_io.c
//  pni_correlation_service
//
//  Created by Benjamin Singer on 10/15/14.
//  Copyright (c) 2014 Benjamin Singer. All rights reserved.
//
#define _GNU_SOURCE
#include <fcntl.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include "pnicorr.h"
#include "pnicorr_fmemopen.h"
#include "pnicorr_io.h"

// =============    helper functions =====================

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
static void save_compressed_result(const float *result, const int M,
                                   const int N, const char *mode,
                                   const int save_text,
                                   const int compression_level, char *outfile) {
  char compress_mode[4];
  sprintf(compress_mode, "%s%d", mode, compression_level);

  if (save_text) {
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
  } else {
    gzFile *fp = gzopen(outfile, compress_mode);
    ERRIF(NULL == fp, "gzopen");
    gzwrite(fp, result, sizeof(float) * M * N);
    gzclose(fp);
  }
}

// ****************** save_result ***************
void pnicorr_savematrix(const float *result, const int M, const int N,
                        const char *mode, const int save_text,
                        const int compression_level, char *outfile) {
  if (compression_level > 0) {
    save_compressed_result(result, M, N, mode, save_text, compression_level,
                           outfile);
  }

  if (save_text) {
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
  } else {
    FILE *fp = fopen(outfile, mode);
    ERRIF(NULL == fp, "fopen");
    fwrite(result, sizeof(float), M * N, fp);
    fclose(fp);
  }
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
  FILE *sp = pnicorr_fmemopen(filebuf, filesize, "r");
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
