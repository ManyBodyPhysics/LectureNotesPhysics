#include "mbpt2V05.hpp"

mbpt2V05::mbpt2V05(twobodyPotential &potential, Model_Space &space) {
  mbpt2V05::name = "mbpt2V05";
  mbpt2V05::potential = &potential;
  mbpt2V05::modelspace = space;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  number_of_elements = pow(number_of_unoccupied_states,2)*pow(number_of_occupied_states,2);
}

double mbpt2V05::getEnergy() {
  double *v = setup_twobody_potential_hhpp();
  double *denominators = setup_vector_of_energy_denominators_hhpp();

  double energy = energy_from_vectors(v, denominators);
  delete v; delete denominators;
  return energy;
}

double * mbpt2V05::setup_twobody_potential_hhpp () {

  double initial_value = 0.0;
  double *v = allocate_and_initialize_double_vector(initial_value);

  int idx = 0;
  int offset = number_of_occupied_states;
  for (int i = 0; i < number_of_occupied_states; ++i) {
    for (int j = 0; j < number_of_occupied_states; ++j) {
      for (int a = 0; a < number_of_unoccupied_states; ++a) {
        for (int b = 0; b < number_of_unoccupied_states; ++b) {
          v[idx++] = 
            potential->get_element(modelspace.qnums, i, j, a+offset, b+offset);
        }
      }
    }
  }
  return v;
}

double * mbpt2V05::allocate_and_initialize_double_vector(double value) {
  double *vector = new double[number_of_elements];
  for ( int i = 0; i < number_of_elements; i++ ) {
    vector[i] = value;
  }
  return vector;
}

double * mbpt2V05::setup_vector_of_energy_denominators_hhpp () {

  double initial_value = 1.0;
  double *denominators = allocate_and_initialize_double_vector(initial_value);


  double *hole_energies = new double[number_of_occupied_states];
  double *particle_energies = new double[number_of_unoccupied_states];

  for (int i = 0; i < number_of_occupied_states; i++) {
    hole_energies[i] = modelspace.qnums[i].energy;
  }

  int offset = number_of_occupied_states;
  for (int i = 0; i < number_of_unoccupied_states; i++) {
    particle_energies[i] = modelspace.qnums[i+offset].energy;
  }

  int idx = 0;
  for (int i = 0; i < number_of_occupied_states; ++i) {
    for (int j = 0; j < number_of_occupied_states; ++j) {
      for (int a = 0; a < number_of_unoccupied_states; ++a) {
        for (int b = 0; b < number_of_unoccupied_states; ++b) {
          denominators[idx++] = hole_energies[i] + hole_energies[j] -
                                particle_energies[a] - particle_energies[b];
        }
      }
    }
  }
  delete hole_energies; delete particle_energies;
  return denominators;
}

double mbpt2V05::energy_from_vectors(double *v, double *denominators) {

  double energy = 0.0;
  for ( int i = 0; i < number_of_elements; i++ ) {
    energy += v[i]*v[i]/denominators[i];
  }
  return 0.25*energy;
}


