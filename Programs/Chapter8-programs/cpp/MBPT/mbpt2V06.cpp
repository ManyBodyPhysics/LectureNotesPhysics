#include "mbpt2V06.hpp"
#include "minnesota_potential.hpp"

mbpt2V06::mbpt2V06(Input_Parameters &parameters, Model_Space &space) {
  mbpt2V06::name = "mbpt2V06";
  mbpt2V06::modelspace = space;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  L = pow((parameters.numberOfParticles)/parameters.density, 1./3.);

}

double mbpt2V06::getEnergy() {
  double energy0;
  double energy = 0.0;
  for(int i = 0; i < modelspace.indhol; ++i){
    for(int j = 0; j < modelspace.indhol; ++j){
      if(i == j){ continue; }
      for(int a = modelspace.indhol; a < modelspace.indtot; ++a){
        for(int b = modelspace.indhol; b < modelspace.indtot; ++b){
          if(a == b){ continue; }
          energy0 = V_Minnesota(modelspace, i, j, a, b, L);
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
