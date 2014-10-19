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

#include "pnicorr.h"

// ******* Memory limits
#define BYTES_PER_GB (1000000000LL)
#define BYTES_PER_MB (1000000LL)

// ******* Helper functions

// ****************** normalize_fvec ***************
static void normalize_fvec(float *data, const int num_timeseries, const int N) {
  const int incX = 1;
  float *dptr = data;
  for (int i = 0; i < num_timeseries; ++i) {
    float norm2 = cblas_snrm2(N, dptr, incX);
    const float alpha = 1.0 / norm2;
    cblas_sscal(N, alpha, dptr, incX);
    dptr += N;
  }
}

// ****************** MAIN ***************
int main(const int argc, const char *argv[]) {

  // ***** parse input args ******
  int normalize = 1;
  pnicorr_iotype_t outtype = pnicorr_iotype_1Dgz;
  long long maxmem = 4000;
  int pre_args = 2;
  const char *ext = "1D.dset";
  const char *comp = ".gz";

  if (argc < pre_args) {
    char *basec = strdup(argv[0]);
    char *bname = basename(basec);
    fprintf(stderr, "\nusage: %s file.1D.dset[.gz] -[no]norm -mem=MB "
                    "-iotype=1D|1Dgz|mat\n",
            bname);
    exit(EXIT_FAILURE);
  }
  const char *filename = argv[1];

  struct opts2struct_t *ops2s = opts2struct_create();
  opts2struct_parseopts(ops2s, argc - pre_args, &argv[pre_args]);

  if (ops2s->found[norm])
    normalize = ops2s->i[norm];
  if (ops2s->found[mem])
    maxmem = ops2s->i[mem];

  if (ops2s->found[iotype]) {
    if (!strncmp(ops2s->iotype, "1D", 2)) {
      outtype = pnicorr_iotype_1D;
      ext = "1D.dset";
      comp = "";
    } else if (!strncmp(ops2s->iotype, "1Dgz", 4)) {
      outtype = pnicorr_iotype_1Dgz;
      ext = "1D.dset";
      comp = ".gz";
    } else if (!strncmp(ops2s->iotype, "mat", 3)) {
      outtype = pnicorr_iotype_mat;
      ext = "mat";
      comp = "";
    }
  }

  //  **********  get data from file ************
  int num_timeseries, trs;
  float *data = pnicorr_load_1D(filename, &num_timeseries, &trs);

  // get outfile name
  const char *outfile = pnicorr_genoutfile(filename, num_timeseries, ext, comp);

  //  **********  Normalize -- unless asked not to via "-nonorm" ************
  if (normalize) {
    LOG("normalizing each row/timeseries independently ...\n");
    TIC;
    normalize_fvec(data, num_timeseries, trs);
    TOC("normalize");
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
  LOG("%lld task(s), since %d^2 floats is %.3f MB and the max MB per "
      "task is %lld MB\n",
      ntasks, num_timeseries, (float)nbytes / (float)BYTES_PER_MB, maxmem);
  float *result = NULL;
  const float alpha = 1.0;
  const float beta = 0.0;

  // get time info immediately
  setbuf(stdout, NULL);

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
    pnicorr_savematrix(result, num_timeseries, num_timeseries, "w", outtype,
                       outfile);
    TOC("saving");
    LOG("done\n");
  } else {
    // ****************  need to do it in blocks and append results *********
    // TODO: use MPI to distribute this, but first let's see how long it takes
    // serially w/OpenMP
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
                         first_task ? "w" : "a", outtype, outfile);
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

  return EXIT_SUCCESS;
}
