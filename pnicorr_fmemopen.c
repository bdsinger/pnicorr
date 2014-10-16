/*

 ported from
 https://github.com/ingenuitas/python-tesseract/blob/master/fmemopen.c

 which is missing; now see
 https://github.com/NimbusKit/memorymapping/blob/master/src/fmemopen.c

 */
#ifndef __linux__

#include <stdlib.h>
#include <string.h>
#include "pnicorr_fmemopen.h"

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

  if (size > available)
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

  if (pos < 0 || pos > mem->size)
    return -1;

  mem->pos = pos;
  return pos;
}

static int closefn(void *handler) {
  free(handler);
  return 0;
}

/* simple, but portable version of fmemopen for OS X / BSD */
FILE *pnicorr_fmemopen(void *buf, size_t size, const char *mode) {
  fmem_t *mem = (fmem_t *)malloc(sizeof(fmem_t));

  memset(mem, 0, sizeof(fmem_t));
  mem->size = size, mem->buffer = buf;
  return funopen(mem, readfn, writefn, seekfn, closefn);
}

/* taviso / rarvmtools */
/*   strchrnul */
/*  for non-linux platforms */
char *pnicorr_strchrnul(const char *s, int c) {
  char *ptr = strchr(s, c);

  if (!ptr) {
    ptr = strchr(s, '\0');
  }

  return ptr;
}

#endif