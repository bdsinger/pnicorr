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
interpolates voxel rows data onto each node (vertex) in a MGH FreeSurfer
surface. But the rows could be any rows data; no spatial information is
used (nor is it supplied by any recognized input arguments).

Normalizes each row by default. Saves correlationss as platform single-precision
floats with .1D.gz output an option. Uses Cblas API so link with the Cblas
implementation library of choice for your platform (Atlas, Gotoblas, MKL,
vecLib, cuBLAS, OpenCL, etc). Tested with GotoBlas2 on 64-bit Linux and vecLib
on 64-bit Mac OS X

Most CBLAS implementations are multicore, so the parallelism (within-machine)
should come from that.

Ben Singer <bdsinger@princeton.edu> April 2013, revised October 2014

*/

#define _GNU_SOURCE
#include <libgen.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <malloc/malloc.h>
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include "pnicorr.h"

// ******* Memory limits
#define BYTES_PER_GB (1000000000LL)
#define BYTES_PER_MB (1000000LL)
#define MAX_FILENAME 256
#define CBLAS_ALPHA (1.0F)
#define CBLAS_BETA (0.0F)
#define PNICORR_ALIGNMENT 64
enum { PNICORR_CORRELATIONS_TAG = 1, PNICORR_NUMTAGS };

#define PNICORR_PROC_BODY_PARAMS                                               \
  (const float *restrict full_matrix, float *restrict submatrix,               \
   float *restrict correlations, const int *restrict offset_for_task,          \
   const int *restrict rows_for_task, const int nt, const int num_rows,        \
   const int num_columns, const char *restrict outfile,                        \
   const pnicorr_iotype_t outtype, const int my_task_id)
#define PNICORR_PROC_BODY_ARGS                                                 \
  (full_matrix, submatrix, correlations, offset_for_task, rows_for_task, nt,   \
   num_rows, num_columns, outfile, outtype, my_task_id)
#define PNICORR_PARAMSET_PARAMS                                                \
  (const int my_task_id, const int *num_rows, const int *num_columns,          \
   long long *ntasks, int *num_workers, float **restrict full_matrix)
#define PNICORR_PARAMSET_ARGS                                                  \
  (my_task_id, &num_rows, &num_columns, &ntasks, &num_workers, &full_matrix)
#define PNICORR_TASK_PROCESSING_PARAMS                                         \
  (const float *restrict full_matrix, float *restrict submatrix,               \
   float *restrict correlations, const int offset_for_task,                    \
   const int rows_for_task, const int num_rows, const int num_columns)
#define PNICORR_TASK_PROCESSING_ARGS(id)                                       \
  (full_matrix, submatrix, correlations, offset_for_task[id],                  \
   rows_for_task[id], num_rows, num_columns)

#define COR_TASK_PROCESSING(id)                                                \
  pnicorr_task_processing PNICORR_TASK_PROCESSING_ARGS(id)

#ifdef HAVE_MPI
/* use MPI */
#include <mpi.h>

void pnicorr_mpi_init(int *restrict my_task_id, int *restrict num_workers);
void pnicorr_mpi_finalize(void);
void pnicorr_mpi_paramset PNICORR_PARAMSET_PARAMS;
void pnicorr_mpi_proc_body PNICORR_PROC_BODY_PARAMS;

#define COR_MPI_INIT(my_task_id, num_workers)                                  \
  pnicorr_mpi_init(&my_task_id, &num_workers)

#define COR_MPI_FINALIZE pnicorr_mpi_finalize()

#define COR_MPI_PARAMSET pnicorr_mpi_paramset PNICORR_PARAMSET_ARGS

#define COR_TASK_PROC_BODY pnicorr_mpi_proc_body PNICORR_PROC_BODY_ARGS

#else
/* no MPI */

#define NOP (void)0
#define COR_MPI_INIT(my_task_id, num_workers) NOP
#define COR_MPI_FINALIZE NOP
#define COR_MPI_PARAMSET NOP
#define COR_TASK_PROC_BODY pnicorr_proc_body PNICORR_PROC_BODY_ARGS

#endif

// ******* Helper functions

// ****************** normalize_fvec ***************
static void normalize_fvec(float *restrict full_matrix, const int num_rows,
                           const int N) {
  const int incX = 1;
  float *dptr = full_matrix;
  for (int i = 0; i < num_rows; ++i) {
    float norm2 = cblas_snrm2(N, dptr, incX);
    const float alpha = 1.0 / norm2;
    cblas_sscal(N, alpha, dptr, incX);
    dptr += N;
  }
}

void pnicorr_task_processing PNICORR_TASK_PROCESSING_PARAMS;

void pnicorr_proc_body PNICORR_PROC_BODY_PARAMS;

void pnicorr_savematrix_for_task(const int i, const int *restrict rows_for_task,
                                 const int num_rows, const char *outfile,
                                 const pnicorr_iotype_t outtype,
                                 const float *restrict correlations);

// ****************** MAIN ***************
int main(int argc, char *argv[]) {
  float *full_matrix = NULL;
  float *correlations = NULL;
  char outfile[MAX_FILENAME];
  int num_rows = 0, num_columns = 0;
  int my_task_id = 0;
  long long ntasks = 1;
  int num_workers = 1;

  COR_MPI_INIT(my_task_id, num_workers);

  // defaults
  int normalize = 1;
  pnicorr_iotype_t outtype = pnicorr_iotype_1Dgz;
  long long maxmem = 4000;

  struct opts2struct_t *ops2s = opts2struct_create();

  if (0 == my_task_id) {

    // ***** parse input args ******
    const char *ext = "1D.dset";
    const char *comp = ".gz";
    const int pre_args = 2;
    char *basec = strdup(argv[0]);
    char *bname = basename(basec);

    if (argc < pre_args) {
      fprintf(stderr, "\nusage: %s file.1D.dset[.gz] -[no]norm -mem=MB "
                      "-iotype=1D|1Dgz|mat\n",
              bname);
      exit(EXIT_FAILURE);
    }
    const char *filename = argv[1];

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

    //  **********  get full_matrix from file ************
    full_matrix = pnicorr_load_1D(filename, &num_rows, &num_columns);

    // get outfile name
    pnicorr_genoutfile(bname, num_rows, ext, comp, outfile);
    free(basec);

    //  **********  Normalize -- unless asked not to via "-nonorm" ************
    if (normalize) {
      LOG("normalizing each row/rows independently ...\n");
      TIC;
      normalize_fvec(full_matrix, num_rows, num_columns);
      TOC("normalize");
    }

    // ************  do the matrix multiplication ***************
    // (results in correlation after normalization)
    long long corbytes =
        (long long)num_rows * (long long)num_rows * (long long)sizeof(float);
    // must split up the task s.t. each task uses MAX_MEM_PER_TASK bytes
    ntasks = (long long)((float)corbytes /
                             (float)((long long)maxmem * BYTES_PER_MB) +
                         0.5);
    LOG("%lld task(s), since %d^2 floats is %.3f MB and the max MB per "
        "task is %lld MB\n",
        ntasks, num_rows, (float)corbytes / (float)BYTES_PER_MB, maxmem);

    if (ntasks <= 1) {
      // get time info immediately
      setbuf(stdout, NULL);

      // *********   all in one go; enough memory to do so *********
      LOG("number of bytes needed for correlations is small enough for a "
          "single "
          "task\n");
      LOG("allocating memory...\n");
      int err =
          posix_memalign((void *)&correlations, PNICORR_ALIGNMENT, corbytes);
      ERRIF(err, "posix_memalign");
      LOG("calling cblas_ssyrk to perform AA' ...\n");

      TIC;
      cblas_ssyrk(CblasRowMajor, CblasUpper, CblasNoTrans, num_rows,
                  num_columns, CBLAS_ALPHA, full_matrix, num_columns,
                  CBLAS_BETA, correlations, num_rows);
      TOC("calculation");

      // *** save ***
      LOG("saving %d x %d task correlations to %s\n", num_rows, num_rows,
          outfile);
      TIC;
      pnicorr_savematrix(correlations, num_rows, num_rows, "w", outtype,
                         outfile);
      TOC("saving");
      LOG("done\n");
    } // ntasks > 1
  }   // my_task_id == 0

  COR_MPI_PARAMSET;

  if (ntasks > 1) {
    float *submatrix = NULL;

    // ****************  need to do it in blocks and append correlations
    // *********
    int nt = (int)ntasks;
    int num_rows_per_task = num_rows / nt;
    if (num_rows_per_task < 1) {
      num_rows_per_task = 1;
      nt = num_rows;
    }

    int *rows_for_task = calloc(num_workers, sizeof(int));
    int *offset_for_task = calloc(num_workers, sizeof(int));
    // have task 0 (master in case of MPI) do the work for leftover rows
    // when total doesn't evenly divide ntasks
    rows_for_task[0] = num_rows - (num_rows_per_task * (nt - 1));
    offset_for_task[0] = 0;
    int maxrows_for_tasks = (rows_for_task[0] > num_rows_per_task)
                                ? rows_for_task[0]
                                : num_rows_per_task;
    for (int i = 1; i < nt; ++i) {
      offset_for_task[i] = offset_for_task[i - 1] + rows_for_task[i - 1];
      rows_for_task[i] = num_rows_per_task;
    }
    int err = posix_memalign((void *)&submatrix, PNICORR_ALIGNMENT,
                             maxrows_for_tasks * num_columns * sizeof(float));
    ERRIF(err, "posix_memalign submatrix");
    err = posix_memalign((void *)&correlations, PNICORR_ALIGNMENT,
                         maxrows_for_tasks * num_rows * sizeof(float));
    ERRIF(err, "posix_memalign correlations");

    COR_TASK_PROC_BODY;

    free(rows_for_task);
    free(offset_for_task);
    if (submatrix) {
      free(submatrix);
      submatrix = NULL;
    }
  }

  if (full_matrix) {
    free(full_matrix);
    full_matrix = NULL;
  }
  if (correlations) {
    free(correlations);
    correlations = NULL;
  }

  COR_MPI_FINALIZE;

  return EXIT_SUCCESS;
}

void pnicorr_task_processing PNICORR_TASK_PROCESSING_PARAMS {

  LOG("submatrix %p full_matrix %p\n", submatrix, full_matrix);

  memcpy(submatrix, full_matrix + offset_for_task,
         rows_for_task * num_columns * sizeof(float));
  memset(correlations, 0, num_rows * rows_for_task * sizeof(float));

  TIC;
  LOG("calling cblas_sgemm using %d num_rows (rows) at a time\n",
      rows_for_task);
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rows_for_task, num_rows,
              num_columns, CBLAS_ALPHA, submatrix, num_columns, full_matrix,
              num_columns, CBLAS_BETA, correlations, num_rows);
  TOC("calculating");
}

void pnicorr_savematrix_for_task(const int i, const int *restrict rows_for_task,
                                 const int num_rows, const char *outfile,
                                 const pnicorr_iotype_t outtype,
                                 const float *restrict correlations) {
  // ******** save/append *********
  LOG("task %d: %s %d x %d task correlations to %s\n", i,
      (i == 0) ? "writing" : "appending", rows_for_task[i], num_rows, outfile);
  TIC;
  pnicorr_savematrix(correlations, rows_for_task[i], num_rows,
                     (i == 0) ? "w" : "a", outtype, outfile);
  TOC("saving");
}

#ifdef HAVE_MPI

void pnicorr_mpi_init(int *my_task_id, int *num_workers) {
  int initialized;
  MPI_Initialized(&initialized);
  if (0 == initialized) {
    MPI_Init(NULL, NULL);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, my_task_id);
  MPI_Comm_size(MPI_COMM_WORLD, num_workers);
}

void pnicorr_mpi_finalize(void) {
  int initialized;
  MPI_Initialized(&initialized);
  if (initialized) {
    int finalized;
    MPI_Finalized(&finalized);
    if (0 == finalized) {
      MPI_Finalize();
    }
  }
}

void pnicorr_mpi_paramset PNICORR_PARAMSET_PARAMS {

  MPI_Bcast((void *)num_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast((void *)num_columns, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast((void *)ntasks, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast((void *)num_workers, 1, MPI_INT, 0, MPI_COMM_WORLD);

  long long fullmatrixbytes = *num_rows * *num_columns * sizeof(float);
  if (my_task_id > 0) {
    int err =
        posix_memalign((void *)full_matrix, PNICORR_ALIGNMENT, fullmatrixbytes);
    ERRIF(err, "posix_memalign");
  } else {
    ALOGIF(*num_workers != *ntasks, "setting ntasks equal to num_workers: (%d "
                                    "workers, %d tasks -> %d tasks)\n",
           *num_workers, (int)*ntasks, *num_workers);
  }
  *ntasks = *num_workers;

  MPI_Bcast((void *)(*full_matrix), *num_rows * *num_columns, MPI_FLOAT, 0,
            MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
}

void pnicorr_mpi_proc_body PNICORR_PROC_BODY_PARAMS {
  COR_TASK_PROCESSING(my_task_id);

  if (0 == my_task_id) {
    pnicorr_savematrix_for_task(0, rows_for_task, num_rows, outfile, outtype,
                                correlations);
    for (int i = 1; i < nt; ++i) {
      MPI_Recv(correlations, num_rows * rows_for_task[i], MPI_FLOAT, i,
               PNICORR_CORRELATIONS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      pnicorr_savematrix_for_task(i, rows_for_task, num_rows, outfile, outtype,
                                  correlations);
    }
  } else {
    LOG("[%d] sending %d x %d correlations\n", my_task_id,
        rows_for_task[my_task_id], num_rows);
    MPI_Send(correlations, num_rows * rows_for_task[my_task_id], MPI_FLOAT, 0,
             PNICORR_CORRELATIONS_TAG, MPI_COMM_WORLD);
  }
}

#else

void pnicorr_proc_body PNICORR_PROC_BODY_PARAMS {
  int num_rows_processed = 0;
  for (int i = 0; i < nt; ++i) {

    COR_TASK_PROCESSING(i);

    pnicorr_savematrix_for_task(i, rows_for_task, num_rows, outfile, outtype,
                                correlations);
    num_rows_processed += rows_for_task[i];

    ALOG("%d%% complete\n",
         (int)((float)num_rows_processed / (float)num_rows * 100.0f));
  }
}

#endif
