/*
        Version 1.0

Correlation service; initial version will accept text input (optionally gzipped)
where comments and any integer (ie containing no ".") in initial columns are
skipped:

[# comments]
[int int int ] fl.oat fl.oat ...

([] means optional.)

Each row is a node, each column a TR, separated by whitespace. This is the
.1D surface fmri data format created by NIH AFNI/SUMA's 3dVol2Surf, which
interpolates voxel timeseries data onto each node (vertex) in a MGH FreeSurfer
surface. But the rows could be any timeseries data; no spatial information is
used (nor is it supplied by any recognized input arguments).

Normalizes each row by default. Saves results as platform single-precision
floats with .1D.gz output an option. Uses Cblas API so link with the Cblas
implementation library of choice for your platform (Atlas, Gotoblas, MKL,
vecLib, cuBLAS, OpenCL, etc). Tested with GotoBlas2 on 64-bit Linux and vecLib
on 64-bit Mac OS X

When GotoBlas2 is compiled with "USE_OPENMP 1" uncommented (in its
Makefile.rules) it will use all available cores in parallel to do the work; no
need to parallelize the for-loop in this code if you use that solution or one
like it. What can    still improve  performance would be the use of multiple
computers (say via MPI) -- for version 2.0

Ben Singer <bdsinger@princeton.edu> April 2013

*/

#define _GNU_SOURCE
#include <libgen.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include "pnicorr_io.h"
#include "pnicorr_reporttime.h"
#include "pnicorr_fmemopen.h"
#include "pnicorr_debugbreak.h"
#include "pnicorr.h"

static const char **flag_text = NULL;
static const char **flag_description = NULL;
static const int *runtime_flags = NULL;

// ******* Helper functions

// ****************** normalize_fvec ***************
static void normalize_fvec(float *X, const int N) {
  const int incX = 1;

  float norm2 = cblas_snrm2(N, X, incX);
  if (norm2 == 0.0 || norm2 == 1.0)
    return;

  const float alpha = 1.0 / norm2;
  cblas_sscal(N, alpha, X, incX);
}

static char get_trailing_arg_char(const char *argv[],
                                  const int *flag_argc_index,
                                  const pnicorr_rflags_t flag) {
  char retfail = (char)-1;

  if (flag_argc_index[flag] == 0)
    return retfail;

  const char *arg = argv[flag_argc_index[flag]];
  size_t prelen = strlen(flag_text[flag]);

  if (strlen(arg) >= prelen + 1) {
    return arg[prelen];
  } else {
    return retfail;
  }
}

// ****************** MAIN ***************
int main(const int argc, const char *argv[]) {
  int compression_level = 0;
  long long maxmem = 4000;

  flag_text = pnicorr_flagtext();
  flag_description = pnicorr_flagdescription();
  runtime_flags = pnicorr_runtimeflags();

  if (argc < 2) {
    char *basec = strdup(argv[0]);
    char *bname = basename(basec);
    fprintf(stderr, "\nusage: %s file.1D[.gz] [options]\n", bname);
    for (pnicorr_rflags_t i = 0; i < pnicorr_rflags_numflags; ++i) {
      fprintf(stderr, " %s:\t%s\n\n", flag_text[i], flag_description[i]);
    }
    exit(EXIT_FAILURE);
  }

  // ****** parse args ********
  const char *filename = argv[1];
  int user_flags = 0;
  int flag_argc_index[pnicorr_rflags_numflags] = {0};
  if (argc > 2) {
    for (int i = 2; i < argc; ++i) {
      for (pnicorr_rflags_t j = 0; j < pnicorr_rflags_numflags; ++j) {
        if (strstr(argv[i], flag_text[j]) != NULL) {
          user_flags |= runtime_flags[j];
          flag_argc_index[j] = i;
          LOG("%s\n", flag_text[j]);
        }
      }
    }
  }

  /*
    user can set the gzip compression level for the output by appending a "1"
    through "9" to the "-gzout" flag as in "-gzout1" through "-gzout9". This
    number is used in the "mode" argument to gzopen, passed to the save_result
    fn
   */
  char trailing_char = (int)get_trailing_arg_char(
      argv, flag_argc_index, pnicorr_rflags_compressoutput);
  compression_level =
      flag_argc_index[pnicorr_rflags_compressoutput] ? DEFAULT_GZIP_LEVEL : 0;
  if ((int)trailing_char > 0) {
    compression_level = (int)trailing_char;
    LOG("compression level set to %d\n", compression_level);
  }

  // -mem variants; determine how much to compute at once
  trailing_char =
      get_trailing_arg_char(argv, flag_argc_index, pnicorr_rflags_memory);
  if ((int)trailing_char >= 0) {
    const char *mm = argv[flag_argc_index[pnicorr_rflags_memory]];
    maxmem = strtol(mm + strlen(flag_text[pnicorr_rflags_memory]), NULL, 10);
    ERRIF(LONG_MAX == maxmem, "strtol");
    LOG("max memory set to %lldM\n", maxmem);
  }

  int num_timeseries, trs;
  //  **********  get data from file ************
  float *data = pnicorr_load_1D(filename, &num_timeseries, &trs);

  //  **********  Normalize -- unless asked not to via "-nonorm" ************
  if (0 == (user_flags & runtime_flags[pnicorr_rflags_nonorm])) {
    LOG("normalizing each row/timeseries independently ...\n");
    float *dptr = data;
    for (int i = 0; i < num_timeseries; ++i) {
      // normalizes each row (timeseries) independently
      normalize_fvec(dptr, trs);
      dptr += trs;
    }
  }

  // ************  do the matrix multiplication ***************
  // (results in correlation after normalization)
  long long int nbytes = (long long int)num_timeseries *
                         (long long int)num_timeseries *
                         (long long int)sizeof(float);
  // must split up the task s.t. each task uses MAX_MEM_PER_TASK bytes
  long long int ntasks =
      (long long int)((float)nbytes /
                          (float)((long long int)maxmem * BYTES_PER_MB) +
                      0.5);
  // LOG("%lld bytes required,%lld MB max, %lld tasks required\n",nbytes,(long
  // long int)maxmem*BYTES_PER_MB,ntasks);
  float *result = NULL;
  const float alpha = 1.0;
  const float beta = 0.0;
  char outfile[MAX_FILENAME];
  char *loc = pnicorr_strchrnul(filename, '.');
  if (*loc != '\0')
    loc = strchr(filename, '.'); // want first one, if it's there at all
  if (loc[1] == '/')
    loc = strchr(loc, '.'); // but if the form ./file.1D, want second one
  strncpy(outfile, filename, loc - filename);
  const char *ext =
      (user_flags & runtime_flags[pnicorr_rflags_savetext]) ? "1D" : "single";
  const char *comp = (compression_level > 0) ? ".gz" : "";
  sprintf(outfile, "%s_%dx%d_correlations.%s%s", outfile, num_timeseries,
          num_timeseries, ext, comp);
  setbuf(stdout, NULL); // get time info immediately

  if (ntasks <= 1) {
    // *********   all in one go; enough memory to do so *********
    LOG("number of bytes needed for result is small enough for a single "
        "task\n");
    LOG("allocating memory...\n");
    result = (float *)calloc((long long int)num_timeseries *
                                 (long long int)num_timeseries,
                             sizeof(float));
    ERRIF(NULL == result, "calloc");
    LOG("calling cblas_ssyrk to perform AA' ...\n");

    TIC;

    cblas_ssyrk(CblasRowMajor, CblasUpper, CblasNoTrans, num_timeseries, trs,
                alpha, data, trs, beta, result, num_timeseries);

    TOC("calculation");

    // *** save ***
    LOG("saving %d x %d task result to %s\n", num_timeseries, num_timeseries,
        outfile);
    TIC;

    pnicorr_savematrix(result, num_timeseries, num_timeseries, "w",
                       user_flags & runtime_flags[pnicorr_rflags_savetext],
                       compression_level, outfile);

    TOC("saving");

    LOG("done\n");
  } else {
    // ****************  need to do it in blocks and append results *********
    // TODO: use MPI to distribute this, but first let's see how long it takes
    // serially w/OpenMP
    LOG("must use %lld tasks, since %d^2 floats is %.3f MB and the max MB per "
        "task is %lld MB\n",
        ntasks, num_timeseries, (float)nbytes / (float)BYTES_PER_MB, maxmem);
    int nt = (int)ntasks;
    int num_timeseries_per_task = num_timeseries / nt;
    if (num_timeseries_per_task < 1) {
      num_timeseries_per_task = 1;
      nt = num_timeseries;
    }
    int leftover_num_timeseries =
        num_timeseries - (num_timeseries_per_task * nt);
    long long int bytes_per_task = (long long int)num_timeseries_per_task *
                                   (long long int)num_timeseries *
                                   (long long int)sizeof(float);
    LOG("allocating %.3f M for task results\n",
        (float)bytes_per_task / (float)BYTES_PER_MB);
    result = (float *)malloc(bytes_per_task);
    ERRIF(NULL == result, "malloc");
    float *offset = data;
    int first_task = 0;
    int last_task = 0;
    int num_timeseries_processed = 0;

    for (int i = 0; i < nt; ++i) {

      first_task = (i == 0);
      last_task = (i == nt - 1);
      if (last_task && (leftover_num_timeseries > 0)) {
        LOG("final task adding %d leftover node remainder\n",
            leftover_num_timeseries);
        num_timeseries_per_task += leftover_num_timeseries;
        bytes_per_task =
            num_timeseries_per_task * num_timeseries * sizeof(float);
        float *new_result = realloc(result, bytes_per_task);
        if (NULL == new_result) {
          free(result);
          ERRIF(1, "realloc");
        }
        result = new_result;
      }

      memset(result, 0, bytes_per_task);

      TIC;

      LOG("calling cblas_sgemm using %d num_timeseries (rows) at a time\n",
          num_timeseries_per_task);
      /*
       multiply this subset by full matrix. ignores the symmetry

       TODO: don't ignore the symmetry. that is: multiply by part of matrix
       with equal or higher numbered columns

       the following tries that, but is not quite right:
       cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans,num_timeseries_per_task,num_timeseries-num_timeseries_processed,trs,alpha,offset,trs,offset,trs,beta,result,num_timeseries);
       */
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                  num_timeseries_per_task, num_timeseries, trs, alpha, offset,
                  trs, data, trs, beta, result, num_timeseries);

      TOC("calculating");

      // ******** save/append *********
      LOG("task %d: %s %d x %d task result to %s\n", i,
          first_task ? "writing" : "appending", num_timeseries_per_task,
          num_timeseries, outfile);
      TIC;

      pnicorr_savematrix(result, num_timeseries_per_task, num_timeseries,
                         first_task ? "w" : "a",
                         user_flags & runtime_flags[pnicorr_rflags_savetext],
                         compression_level, outfile);

      TOC("saving");

      if (!last_task) {
        offset += (trs * num_timeseries_per_task);
      }
      num_timeseries_processed += num_timeseries_per_task;

      ALOG("%d%% complete\n", (int)((float)num_timeseries_processed /
                                    (float)num_timeseries * 100.0f));
    }
    assert(num_timeseries_processed == num_timeseries);
  }

  if (data)
    free(data);
  if (result)
    free(result);
  if (runtime_flags)
    free((void *)runtime_flags);
  if (flag_text)
    free(flag_text);
  if (flag_description)
    free(flag_description);

  return EXIT_SUCCESS;
}
