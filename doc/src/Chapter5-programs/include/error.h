#include <string>
#include <stdlib.h>
using namespace std;

// Header file for the verbose class. An object of this class
// called ERR should be created at the highest scope (outside 
// main). The header file declares ERR as external.

#ifndef INCLUDED_ERROR
#define INCLUDED_ERROR
class Error
{
  public:
    void FileR(const char*, const char*);
    void FileW(const char*, const char*);
    void FileA(const char*, const char*);
    void General(const char*, const char*, ...);
    void NotImplemented(const char*, const char*, ...);

};

extern Error ERR;
#endif
