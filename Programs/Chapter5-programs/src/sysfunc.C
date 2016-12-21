#include <sysfunc.h>
#include <cstddef>
#include <stdlib.h>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

//---- Some file scoped variables
namespace Comms
{
  static int rank = -1;
  static int size = -1;
}

bool Comms::initialized = false;

void Comms::Initialize()
{

  if (initialized) { return; }

#ifdef USE_MPI

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#else

  rank = 0;
  size = 1;

#endif

  initialized = true;

}

void Comms::Finalize() {

  if (!initialized) { return; }

#ifdef USE_MPI
  MPI_Finalize();
#else
#endif

}


void Comms::RaiseError() {
  Finalize();
  exit(EXIT_FAILURE);
}

int Comms::Rank()
{
  if (!initialized) { Initialize(); }
  return rank;
}

int Comms::Size()
{
  if (!initialized) { Initialize(); }
  return size;
}


void Comms::Sync() {
#ifdef USE_MPI
  if (!Comms::initialized) { Comms::Initialize(); }
  MPI_Barrier(MPI_COMM_WORLD);
#else
#endif
}



complex<double> Comms::GlobalSum(complex<double> z)
{

  if (!initialized) { Initialize(); }

#ifdef USE_MPI
  double sendbuf[2];
  double recvbuf[2];

  sendbuf[0] = z.real();
  sendbuf[1] = z.imag();

  MPI_Reduce( sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

  return complex<double>(recvbuf[0], recvbuf[1]);
#else
  return z;
#endif


}
