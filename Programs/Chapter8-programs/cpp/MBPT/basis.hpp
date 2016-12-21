#ifndef BASIS_H
#define BASIS_H
  typedef struct basis basis;

  struct basis {
    int *t;
    int *m;
    int *nx;
    int *ny;
    int *nz;
    double *energy;
  };
#endif
