#include <string>
#include <iostream>
#include <stdarg.h>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
using namespace std;
#include "error.h"
#include "sysfunc.h"

Error ERR;

void Error::FileR(const char* fname, const char* filename)
{

   printf("%s : Cannot open file %s to read.\n", fname, filename);
  Comms::RaiseError();

}

void Error::FileW(const char* fname, const char* filename)
{

    printf("%s : Cannot open file %s to write.\n", fname, filename);
  Comms::RaiseError();

}

void Error::FileA(const char* fname, const char* filename)
{


    printf("%s : Cannot open file %s to append.\n", fname, filename);
  Comms::RaiseError();

}

void Error::General(const char* fname, const char* format, ...)
{

  va_list args;
  va_start(args, format);
    printf("%s : General error : ", fname);
    vprintf(format, args);
    printf("\n");
  Comms::RaiseError();

}

void Error::NotImplemented(const char* fname, const char* format, ...)
{

  va_list args;
  va_start(args, format);
    printf("%s : Not implemented : ", fname);
   vprintf(format, args);
   printf("\n");
  Comms::RaiseError();

}
