#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <time.h>
using namespace std;
#include "constants.h"
#include "enum.h"
#include "arg.h"
#include "lattice.h"
#include "propagator.h"
#include "random.h"
#include "one_body.h"
#include "two_body.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  char* fname = "int main()";
  char* prog_name = "two_body";

  VRB.Flow(fname, "Running production code %s, compiled on %s at %s. ", prog_name, __DATE__, __TIME__ );
  VRB.LibInfo(fname);

  //---- Establish communications
  Comms::Initialize();
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg;
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg); 

  LatticeArg lattice2_arg;
  lattice2_arg.Decode("args/lattice2.arg");

  EvoArg evo_arg;
  evo_arg.Decode("args/evo.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  InteractionArg interaction_arg;
  interaction_arg.Decode("args/interaction.arg");

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  PropagatorArg prop_arg;
  prop_arg.Decode("args/prop.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  TwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");

  //---- Perform some checks
  if (evo_arg.unload_period<1) { ERR.General(fname,"evo_arg.unload_period must be greater than zero."); }

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattice (contains two and three-body fields, initialized when instantiated)
  Lattice lattice2(lattice2_arg, &rng);

  //---- Instantiate two- and three-body hamiltonians, as well as the external potential
  Interaction3 interaction(&lattice2, interaction_arg, kinetic_arg, kinetic_arg, 1.0);
  Potential potential(potential_arg, 0.5);
  Kinetic kinetic(kinetic_arg, 0.5);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

  //---- Some measurement classes for evaluating two and three-body correlators (assumes different species)

  OneBody source1(one_body_arg);
  OneBody source2(one_body_arg);
  source1.Set("args/source1.arg");
  source2.Set("args/source2.arg");
  TwoBody two_body(two_body_arg);

  //---- Open files for measurements
  //---- row = config number, col = time separation
  FILE* file;
  char f_name[99];
  sprintf(f_name, "%s", (two_body_arg.file_stem).c_str()); 
  file = Fopen(IO_TYPE_ROOT,f_name,"a");

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<double> corr[GJP.Tsites()];
  for (int t=0; t<GJP.Tsites(); t++) { corr[t]=0.0; }

  //---- Instantiate propagator classes
  Propagator prop1(prop_arg);
  Propagator prop2(prop_arg);

  clock_t clock_base = clock();
  double  clock_elapsed;

  //---- Begin the simulation by looping over configurations
  for (int j=0; j<configurations; j++) {

    clock_elapsed = double(clock()-clock_base)/CLOCKS_PER_SEC;
    VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);
    VRB.Flow(fname, "Configuration %d of %d (with an unload_period of %d).",j+1,configurations, unload_period);

    //---- Save the current RNG state
    rng.WriteState(random_arg.file_stem.c_str());

    for (int k=0; k<unload_period; k++) {

      //---- Set the propagator sources and inititialize
      prop1.SetSource(source1.Get());
      prop2.SetSource(source2.Get());
    
      for(int t=0; t<GJP.Tsites(); t++) {

//        VRB.Debug(fname, "t = %d", t);
    
        //---- Perform contractions
        corr[t] += two_body.Run( prop1.Get(), prop2.Get() );

        //---- Compute propagators on the current background configuration
        prop1.Run(hamiltonians);
        prop2.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice2.Refresh();

      }
    }

    //---- Store accumulated results in file 
    complex<double> corr_glb_sum;
    for(int t=0; t<GJP.Tsites(); t++) {

      //---- Need to perform a global average across all nodes
      corr_glb_sum = Comms::GlobalSum(corr[t]/(double)unload_period);
      corr_glb_sum /= Comms::Size();

      //---- Then write the result to a file
      Fprintf(file, "%.*e %.*e ", PREC, corr_glb_sum.real(), PREC, corr_glb_sum.imag());

      //---- Reset accumulated corr
      corr[t] = 0.0;
    }

    Fprintf(file, "\n");
    Fflush(file);

  }

  //---- Close data file
  Fclose(file);

  Comms::Finalize();

  return(EXIT_SUCCESS);
}
