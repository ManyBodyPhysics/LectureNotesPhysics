#include <string>
#include "arg.h"
using namespace std;

// Header file for the verbose class. An object of this class
// called VRB should be created at the highest scope (outside 
// main). The header file declares VRB as external.

#ifndef INCLUDED_VERBOSE
#define INCLUDED_VERBOSE
class Verbose
{
  private:
    bool func_level; // function level switch
    bool warn_level; // warning level switch
    bool result_level; // result level switch
    bool flow_level; // flow level switch
    bool debug_level; // debug level switch

  public:
    Verbose();
    ~Verbose();

    void SetLevel(VerboseArg);
    void SetFuncLevel(bool);
    void SetWarnLevel(bool);
    void SetResultLevel(bool);
    void SetFlowLevel(bool);
    void SetDebugLevel(bool);

    void Func(const char*);
    void Warn(const char*, const char*, ...);
    void Result(const char*, const char*, ...);
    void Flow(const char*, const char*, ...);
    void Debug(const char*, const char*, ...);

    void LibInfo(const char* fname);

};

extern Verbose VRB;
#endif
