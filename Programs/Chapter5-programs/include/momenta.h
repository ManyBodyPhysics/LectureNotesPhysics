#include <vector>
using namespace std;
#include "enum.h"
#include "arg.h"

#ifndef INCLUDED_MOMENTA
#define INCLUDED_MOMENTA

class Momenta
{
  private:
    double* omega;
    vector<int>* momentum;
    vector<int> shell_count;
    int* parity_list;
    int momenta_count;
  public:
    Momenta(MomentaArg);
    ~Momenta();
    int MomentaCount();
    vector<int> GetMomentum(int);
    int OppositeParityIndex(int);
    int ShellCount(int);
    int AccShellCount(int);
    int NumShells();
};

#endif
