#include "enum.h"
#include <string>
#include <vector>
using namespace std;
#ifndef INCLUDED_ARG
#define INCLUDED_ARG

class DoArg
{
  public:
    DoArg();
    void Decode(string); // class function that imports parameters from a file
    int t_sites; // number of sites in the t direction
    int x_sites; // number of sites in the x direction
    int y_sites; // number of sites in the y direction
    int z_sites; // number of sites in the z direction
    BoundaryType x_boundary_type; // boundary conditions in the x direction for the fermion
    BoundaryType y_boundary_type; // boundary conditions in the y direction for the fermion
    BoundaryType z_boundary_type; // boundary conditions in the z direction for the fermion
    double cutoff;
};

class LatticeArg
{
  public:
    LatticeArg();
    void Decode(string); // class function that imports parameters from a file
    enum FieldType field_type; // type of two-body auxiliary field
};

class RandomArg
{
  public:
    RandomArg();
    void Decode(string); // class function that imports parameters from a file
    enum RandomType random_type; // type of random number generator
    enum SeedType seed_type; // type of seed for the random number generator
    long seed; // input value for the seed
    string file_stem; // file stem for random state 
};

class PropagatorArg
{
  public:
    PropagatorArg();
    void Decode(string); // class function that imports parameters from a file
    int n_species; // number of species/
    string file_stem; // file name stem
};

class OneBodyArg
{
  public:
    OneBodyArg();
    void Decode(string); // class function that imports parameters from a file
    SourceType source_type; // one-particle basis state type for the source
    double lambda1;
    double lambda2;
    double lambda3;
};

class TwoBodyArg
{
  public:
    TwoBodyArg();
    void Decode(string); // class function that imports parameters from a file
    WavefuncType wavefunc_type; // two-body wave function type
    DispersionType dispersion_type1;
    DispersionType dispersion_type2;
    double mass1; // wave function mass parameter for first particle
    double mass2; // wave function mass parameter for second particle
    double lambda; // additional wave function parameter 
    string file_stem; // file name stem
};

class VerboseArg
{
  public:
    VerboseArg();
    void Decode(string);
    bool func_level; // function level switch
    bool warn_level; // warning level switch
    bool result_level; // result level switch
    bool flow_level; // flow level switch
    bool debug_level; // debug level switch
};

class EvoArg 
{
  public:
    EvoArg();
    void Decode(string);
    int start; // configuration start number
    int unload_period; // configuration unload period
    int configurations; // total number of configurations to be generated
};

class MomentaArg
{
  public:
    MomentaArg();
    void Decode(string);
    DispersionType dispersion_type; // dispersion relation used to compute fermi-energy
    double fermi_energy; // fermi energy
    double mass; // fermion mass used to compute momenta below fermi energy
};

class InteractionArg
{
  public:
    InteractionArg();
    void Decode(string);
    InteractionType interaction_type;
    vector<double> couplings;
    int num_couplings;
};

class PotentialArg
{
  public:
    PotentialArg();
    void Decode(string);
    PotentialForm potential_form; // Potential functional form (none, harmonic, coulomb)
    PotentialType potential_type; // potential discretization type
    double spring_constant1; // harmonic potential spring constant in x-direction
    double spring_constant2; // harmonic potential spring constant in y-direction
    double spring_constant3; // harmonic potential spring constant in z-direction
};

class KineticArg
{
  public:
    KineticArg();
    void Decode(string); // class function that imports parameters from a file
    DispersionType dispersion_type; // Kinetic operator dispersion relation
    double mass; // mass for D operator
};


#endif
