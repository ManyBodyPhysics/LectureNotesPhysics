#include "mbpt2V00.hpp"

mbpt2V00::mbpt2V00(twobodyPotential &potential, Model_Space &space) {
  mbpt2V00::name = "mbpt2V00";
  mbpt2V00::potential = &potential;
  mbpt2V00::modelspace = space;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
}

double mbpt2V00::getEnergy() {
  double energy0;
  double energy = 0.0;
  for(int i = 0; i < modelspace.indhol; ++i){
    for(int j = 0; j < modelspace.indhol; ++j){
      if(i == j){ continue; }
      for(int a = modelspace.indhol; a < modelspace.indtot; ++a){
        for(int b = modelspace.indhol; b < modelspace.indtot; ++b){
          if(a == b){ continue; }
          energy0 = potential->get_element(modelspace.qnums, i, j, a, b);
          energy0 *= energy0;
          energy0 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                          modelspace.qnums[a].energy - modelspace.qnums[b].energy);
          energy += energy0;
        }
      }
    }
  }
  energy *= 0.25;
  return energy;
}
