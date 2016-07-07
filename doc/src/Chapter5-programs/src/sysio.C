#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <cstring>
#include "sysio.h"
#include "sysfunc.h"
#include "error.h"

static FILE FAKE;
static FILE *FAKE_P = &FAKE;
const int MAX_FILENAME =200;

FILE *Fopen(IoType type, const char *filename, const char *mode) {

  char* fname = "FILE *Fopen(IoType, const char*, const char *)";

  FILE *fp = NULL;
  if ( type == IO_TYPE_ROOT && Comms::Rank() ) return &FAKE;
  if(type == IO_TYPE_ALL){
    char f_name[MAX_FILENAME];
    if(strlen(filename)+6 > MAX_FILENAME) {
      ERR.General(fname,"File name (%s) is too long.",filename);
      return NULL;
    }
    sprintf(f_name,"%s.%d", filename, Comms::Rank());
    fp = fopen(f_name,mode);
  } else {
    fp =  fopen(filename,mode);
  }
  if (fp==NULL)
    ERR.General(fname,"Cannot open %s",filename);
  return fp;
}

int Fclose(FILE *stream) {
  if ( stream == &FAKE )  return 1;
  return fclose(stream);
}

int Fprintf(FILE *stream, const char *format,...) {
  if ( stream == &FAKE )  return 1;
  va_list args;
  va_start(args,format);
  int nb = vfprintf(stream, format, args);
  va_end(args);
  return nb;
}

size_t Fwrite(const void *ptr, size_t size, size_t n, FILE *stream){
  if ( stream == FAKE_P )  return n;
  return fwrite(ptr, size, n, stream);
}

size_t Fread(void *ptr, size_t size, size_t n, FILE *stream){
  char* fname = "size_t Fread(void *, size_t, size_t, FILE*)";
  if ( stream == FAKE_P ) {
    ERR.General(fname,"Trying to read from invalid file stream, possibly from using Fopen with wrong flags");
  }
  return fread(ptr, size, n, stream);
}

int Fflush(FILE *stream) {
  if ( stream == &FAKE )  return NULL;
  return fflush(stream);
}


