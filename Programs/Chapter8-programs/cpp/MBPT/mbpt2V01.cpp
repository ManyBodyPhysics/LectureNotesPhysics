#include "mbpt2V01.hpp"

mbpt2V01::mbpt2V01(twobodyPotential &potential, Model_Space &space) {
  mbpt2V01::name = "mbpt2V01";
  mbpt2V01::potential = &potential;
  mbpt2V01::modelspace = space;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
}

double mbpt2V01::getEnergy() {
  double energy = 0.0;
  #pragma omp parallel shared(energy)
  {
    double temp;
    #pragma omp for reduction(+:energy)
    for(int i = 0; i < modelspace.indhol; ++i){
      for(int j = 0; j < modelspace.indhol; ++j){
        if(i == j){ continue; }
        for(int a = modelspace.indhol; a < modelspace.indtot; ++a){
          for(int b = modelspace.indhol; b < modelspace.indtot; ++b){
            if(a == b){ continue; }
            temp = potential->get_element(modelspace.qnums, i, j, a, b);
            temp *= temp;
            temp /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                            modelspace.qnums[a].energy - modelspace.qnums[b].energy);
            energy += temp;
          }
        }
      }
    }
  }
  return 0.25*energy;
}
