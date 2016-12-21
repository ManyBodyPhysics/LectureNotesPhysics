#include <fourier.h>
#include "global_job_parameter.h"
#include "verbose.h"
#include "error.h"

//---- Some file scoped variables
namespace Fourier
{
  static fftw_plan p1;
  static fftw_plan p2;
}

int Fourier::instances = 0;
fftw_complex* Fourier::in;
fftw_complex* Fourier::out;
Phase Fourier::phase;

void Fourier::Initialize()
{

  char* fname = "void Fourier::Initialize()";

  if (instances==0) {
    //---- This is the first instance of Fourier
    VRB.Debug(fname, "Allocating memory and creating plans for FFTW.");
    in = (fftw_complex*) fftw_malloc(GJP.Vol()*sizeof(fftw_complex));
    out = (fftw_complex*) fftw_malloc(GJP.Vol()*sizeof(fftw_complex));
    p1 = fftw_plan_dft_3d(GJP.Xsites(), GJP.Ysites(), GJP.Zsites(), in, out, FFTW_BACKWARD, FFTW_EXHAUSTIVE); 
    p2 = fftw_plan_dft_3d(GJP.Xsites(), GJP.Ysites(), GJP.Zsites(), in, out, FFTW_FORWARD, FFTW_EXHAUSTIVE);
  }

  instances++;

}

void Fourier::Finalize()
{

  char* fname = "void Fourier::Finalize()";

  instances--;
  
  if (instances==0) {
    //---- This is the last instance of Hamiltonian
    VRB.Debug(fname, "Deallocating memory and destoying plans for FFTW.");
    fftw_free(in);
    fftw_free(out);
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_cleanup();
  }

}

void Fourier::Forward()
{
  fftw_execute(p2);
}

void Fourier::Backward()
{
  fftw_execute(p1);
}
