#include <iomanip>
#include "mbpt_factory.hpp"
#include "mbSolver.hpp"
#include "timer.hpp"

int main (int argc, char * argv[]) {

  std::cout << std::setprecision(12);

  Input_Parameters Parameters = parse_commandline_arguments(argc, argv);

  mbSolver *solver = mbptFactory(Parameters);
  if (!solver) {
    std::cout << "Not a valid mbpt type!" << std::endl;
    return 1;
  }
  std::cout << "Selected method: " << solver->getName() << std::endl;

  timer time = timer();
  time.start();
  double energy = solver->getEnergy();
  double elapsed = time.get_elapsed();

  std::cout << "Energy: " << energy << std::endl;
  std::cout << "Time to run method: " << elapsed << std::endl;
  
  return 0;
}
