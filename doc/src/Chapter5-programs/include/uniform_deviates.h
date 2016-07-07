#include "arg.h"

#ifndef INCLUDED_UNIFORM_DEVIATES
#define INCLUDED_UNIFORM_DEVIATES
class UniformDeviate
// Abstract random number generator base class
{
  public:
    UniformDeviate() {};
    virtual ~UniformDeviate() {};
    virtual void Init(long)=0; // Initialize the generator
    virtual double Run()=0; // Generate a random number
    virtual int StateSize()=0; // Determine the size of the generator state
    virtual void GetState(long*)=0; // Save the state in integer array
    virtual void SetState(long*)=0; // Set the state from integer array
};


class Ran0 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    long idum;

  public:
    Ran0();
    ~Ran0();
    void Init(long);
    double Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
};

class Ran1 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    long idum;
    long iy;
    long* iv;

  public:
    Ran1();
    ~Ran1();
    void Init(long);
    double Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
};



class Ran2 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    long idum;
    long idum2;
    long iy;
    long* iv;

  public:
    Ran2();
    ~Ran2();
    void Init(long);
    double Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
};

class Ran3 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    int inext;
    int inextp;
    long* ma;

  public:
    Ran3();
    ~Ran3();
    void Init(long);
    double Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
};

//---- Ranlxs random generator declarations
typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

class Ranlxs : public UniformDeviate
// Random number generator class derived from Random
{
  private:
//    int init;
    int pr;
    int prm;
    int ir;
    int jr;
    int is;
    int is_old;
    int next[96];
    float one_bit;
    vec_t carry;
    union { dble_vec_t vec[12]; int num[96]; } x;
 
    void Update(void);
    void DefineConstants(void);

    void LuxLevel(int);

  public:
    Ranlxs(int);
    ~Ranlxs();
    void Init(long);
    double Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
};

class Ranlxd : public UniformDeviate
// Random number generator class derived from Random
{
  private:
//    int init;
    int pr;
    int prm;
    int ir;
    int jr;
    int is;
    int is_old;
    int next[96];
    float one_bit;
    vec_t carry;
    union { dble_vec_t vec[12]; int num[96]; } x;
 
    void Update(void);
    void DefineConstants(void);

    void LuxLevel(int);

  public:
    Ranlxd(int);
    ~Ranlxd();
    void Init(long);
    double Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
};

#endif
