#include "mbpt3V00.hpp"
#include "minnesota_potential.hpp"

mbpt3V00::mbpt3V00(Input_Parameters &parameters, Model_Space &space) {
  mbpt3V00::name = "mbpt3V00";
  mbpt3V00::modelspace = space;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  number_of_states = number_of_occupied_states + number_of_unoccupied_states;
  L = pow((parameters.numberOfParticles)/parameters.density, 1./3.);

}

double mbpt3V00::getEnergy() {
  double energy0, energy1;
  double energy = 0.0;
  for(int i = 0; i < number_of_occupied_states; ++i){
    for(int j = 0; j < number_of_occupied_states; ++j){
      if(i == j){ continue; }
      for(int a = number_of_occupied_states; a < number_of_states; ++a){
        for(int b = number_of_occupied_states; b < number_of_states; ++b){
          if(a == b){ continue; }
          energy0 = V_Minnesota(modelspace, i, j, a, b, L);
          energy0 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                            modelspace.qnums[a].energy - modelspace.qnums[b].energy);
          for(int c = number_of_occupied_states; c < number_of_states; ++c){
            for(int d = number_of_occupied_states; d < number_of_states; ++d){
              if(c == d){ continue; }
              energy1 = V_Minnesota(modelspace, a, b, c, d, L);
              energy1 *= V_Minnesota(modelspace, c, d, i, j, L);
              energy1 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                            modelspace.qnums[c].energy - modelspace.qnums[d].energy);
              energy += energy0 * energy1;
            }
          }
        }
      }
    }
  }
  energy *= 0.125;

  return energy;
}
