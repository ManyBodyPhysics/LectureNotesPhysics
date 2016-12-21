#include <time.h>
#include "arg.h"
#include "random.h"
#include "uniform_deviates.h"
#include <math.h>
#include <error.h>
#include <verbose.h>
#include <sysfunc.h>
#include <sysio.h>

Random::Random(RandomArg random_arg)
{
  char* fname = "Random::Random(RandomArg)";

  random_type = random_arg.random_type;

  switch(random_arg.random_type) 
  {
    case RANDOM_TYPE_RAN0:
      VRB.Flow(fname,"Using RAN0");
      r =  new Ran0;
      break;
    case RANDOM_TYPE_RAN1:
      VRB.Flow(fname,"Using RAN1");
      r =  new Ran1;
      break;
    case RANDOM_TYPE_RAN2:
      VRB.Flow(fname,"Using RAN2");
      r =  new Ran2;
      break;
    case RANDOM_TYPE_RAN3:
      VRB.Flow(fname,"Using RAN3");
      r =  new Ran3;
      break;
    case RANDOM_TYPE_RANLXS0:
      VRB.Flow(fname,"Using RANLXS0");
      r =  new Ranlxs(0);
      break;
    case RANDOM_TYPE_RANLXS1:
      VRB.Flow(fname,"Using RANLXS1");
      r =  new Ranlxs(1);
      break;
    case RANDOM_TYPE_RANLXS2:
      VRB.Flow(fname,"Using RANLXS2");
      r =  new Ranlxs(2);
      break;
    case RANDOM_TYPE_RANLXD1:
      VRB.Flow(fname,"Using RANLXD1");
      r =  new Ranlxd(1);
      break;
    case RANDOM_TYPE_RANLXD2:
      VRB.Flow(fname,"Using RANLXD2");
      r =  new Ranlxd(2);
      break;
    default:
      ERR.NotImplemented(fname, "random_type.");
  }

  //---- Initialize the random number generator
  switch(random_arg.seed_type) 
  {
    case SEED_TYPE_DEFAULT:
      VRB.Flow(fname,"Using DEFAULT as seed.");
      random_arg.seed = 1;
      break;  
    case SEED_TYPE_INPUT:
      VRB.Flow(fname,"Using INPUT as seed.");
      if (random_arg.seed <= 0) ERR.General(fname, "Input seed must be positive.");
      break;  
    case SEED_TYPE_TIME:
      VRB.Flow(fname,"Using TIME as seed.");
      random_arg.seed = time((time_t *)NULL);
      break;  
    case SEED_TYPE_FILE:
      VRB.Flow(fname,"Using FILE as seed.");
      break;  
    default:
      ERR.NotImplemented(fname, "Unrecognized seed_type.");
  }

  if (random_arg.seed_type==SEED_TYPE_FILE) {
    ReadState(random_arg.file_stem.c_str());
  } else {
    random_arg.seed += 491*Comms::Rank();
    VRB.Flow(fname,"Initializing RNG seed (%d/%d) is %ld\n", Comms::Rank(), Comms::Size()-1, random_arg.seed);
    r->Init(random_arg.seed);
  }

}

Random::~Random()
{
  delete r;
}

double Random::Uniform()
{
  return r->Run();
}

double Random::Gauss(double sigma)
{
  double x, y, r2;

  do {
    x = 2*Uniform() - 1;
    y = 2*Uniform() - 1;;
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0);

  // Box-Muller transform
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

int Random::Z(int p)
{
  return (int)floor(p*Uniform());
}

void Random::WriteState(const char* f_name)
{
  char* fname = "Random::WriteState(const char*)";
  VRB.Func(fname);

  long *state;
  int state_size = r->StateSize();
  state = (long *)malloc(state_size*sizeof(long));
  size_t count;
  FILE* file;

  VRB.Flow(fname, "Writing state to file: %s",f_name);
  r->GetState(state);
  file = Fopen(IO_TYPE_ALL, f_name, "wb");
  count = Fwrite(&random_type, sizeof(RandomType), 1, file);
  if ((int)count!=1) { ERR.General(fname, "Cannot write sead type to file."); }
  count = Fwrite(state, sizeof(long), state_size, file);
  if ((int)count!=state_size) { ERR.General(fname, "Cannot write state to file."); }
  Fclose(file);

  free(state);
}

void Random::ReadState(const char* f_name)
{
  char* fname = "Random::ReadState(const char*)";
  VRB.Func(fname);

  long *state;
  int state_size = r->StateSize();
  state = (long *)malloc(state_size*sizeof(long));
  size_t count;
  FILE* file;
  RandomType file_random_type;

  VRB.Flow(fname, "Reading state from file: %s",f_name);
  file = Fopen(IO_TYPE_ALL, f_name,"rb");
  count = Fread(&file_random_type, sizeof(RandomType), 1, file);
  if ((int)count!=1) { ERR.General(fname, "Cannot read sead type."); }
  if (random_type!=file_random_type) { ERR.General(fname, "Random type of file does not match."); }
  count = Fread(state, sizeof(long), state_size, file);
  if ((int)count!=state_size) { ERR.General(fname, "File may be corrupted."); }
  Fclose(file);
  r->SetState(state);

  free(state);

}


