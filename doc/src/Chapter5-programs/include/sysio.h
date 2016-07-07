#include "enum.h"
#include "config.h"

#ifndef INCLUDED_SYSIO
#define INCLUDED_SYSIO


FILE *Fopen(IoType type, const char *filename, const char *mode);
int Fclose(FILE *stream);
int Fprintf(FILE *stream, const char *format, ...);
size_t Fwrite(const void *ptr, size_t size, size_t n, FILE *stream);
size_t Fread(void *ptr, size_t size, size_t n, FILE *stream);
int Fflush(FILE *stream);

#endif
