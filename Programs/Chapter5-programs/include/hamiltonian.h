#include "arg.h"
#include "lattice.h"

#ifndef INCLUDED_HAMILTONIAN
#define INCLUDED_HAMILTONIAN

// Abstract class from which HamiltonianYukawa and Potential are derived
class Hamiltonian
{
  private:

  public:
    Hamiltonian();
    virtual ~Hamiltonian();
    virtual void Evolve(double*) = 0;
};

// Classes derived from Hamiltonian 

class Interactions : public Hamiltonian
{
  private:
    int vol;
    double *a;
    double *c;
    vector<Hamiltonian*> interactions;
  public:
    Interactions( vector<Hamiltonian*> );
    ~Interactions();
    void Evolve(double*);
};

class Interaction : public Hamiltonian
{
  private:
    Lattice* lattice;
    double* interaction;
  public:
    Interaction(Lattice*, InteractionArg, KineticArg, KineticArg, double);
    ~Interaction();
    void Evolve(double*);
};

class Interaction2 : public Hamiltonian
{
  private:
    Lattice* lattice;
    double* tmp;
    double* interaction;
  public:
    Interaction2(Lattice*, InteractionArg, KineticArg, KineticArg, double);
    ~Interaction2();
    void Evolve(double*);
};

class Interaction3 : public Hamiltonian
{
  private:
    Lattice* lattice;
    double* tmp;
    double* interaction;
  public:
    Interaction3(Lattice*, InteractionArg, KineticArg, KineticArg, double);
    ~Interaction3();
    void Evolve(double*);
};

class Potential : public Hamiltonian
{
  private:
    double* potential;
    PotentialForm potential_form;
    PotentialType potential_type;
  public:
    Potential(PotentialArg, double);
    ~Potential();
    void Evolve(double*);
};

class Kinetic : public Hamiltonian
{
  private:
    double* xi;
  public:
    Kinetic(KineticArg, double);
    ~Kinetic();
    void Evolve(double*);
};



#endif
