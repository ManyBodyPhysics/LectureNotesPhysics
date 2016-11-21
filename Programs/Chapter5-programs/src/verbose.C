#include <string>
#include <iostream>
#include <stdarg.h>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
using namespace std;
#include "verbose.h"
#include "arg.h"
#include "config.h"

Verbose VRB;

Verbose::Verbose()
{

  char* fname = "Verbose::Verbose()";

  func_level = true;
  warn_level = true;
  result_level = true;
  flow_level = true;
  debug_level = true;

}

Verbose::~Verbose()
{
}


void Verbose::SetFuncLevel(bool b) { func_level = b; }
void Verbose::SetWarnLevel(bool b) { warn_level = b; }
void Verbose::SetResultLevel(bool b) { result_level = b; }
void Verbose::SetFlowLevel(bool b) { flow_level = b; }
void Verbose::SetDebugLevel(bool b) { debug_level = b; }

void Verbose::SetLevel(VerboseArg verbose_arg)
{
  SetFuncLevel(verbose_arg.func_level);
  SetWarnLevel(verbose_arg.warn_level);
  SetResultLevel(verbose_arg.result_level);
  SetFlowLevel(verbose_arg.flow_level);
  SetDebugLevel(verbose_arg.debug_level);
}

void Verbose::Func(const char* fname)
{
  if (func_level==true) {
     printf("%s : Entered :\n", fname);
  }
}

void Verbose::Warn(const char* fname, const char* format, ...)
{

  if (warn_level==true) {
    va_list args;
    va_start(args, format);
     printf("%s : Warning : ", fname);
     vprintf(format, args);
     printf("\n");
  }
}

void Verbose::Result(const char* fname, const char* format, ...)
{

  if (result_level==true) {
    va_list args;
    va_start(args, format);
      printf("%s : Result : ", fname);
      vprintf(format, args);
      printf("\n");
  }
}

void Verbose::Flow(const char* fname, const char* format, ...)
{

  if (flow_level==true) {
    va_list args;
    va_start(args, format);
     printf("%s : Flow : ", fname);
      vprintf(format, args);
     printf("\n");
  }
}


void Verbose::Debug(const char* fname, const char* format, ...)
{

  if (debug_level==true) {
    va_list args;
    va_start(args, format);
    // printf("%s : Debug : ", fname);
    // vprintf(format, args);
    // printf("\n");
  }
}

void Verbose::LibInfo(const char* fname)
{
  //printf("%s : Verbose::LibInfo(const char*) : Using %s, compiled on %s at %s.\n", fname, VERSION_STR,__DATE__, __TIME__);
}
