#include "MBPTfunctions.hpp"

int main(int argc, char * argv[])
{
  Input_Parameters Parameters;
  Model_Space Space;
  Channels Chan;

  // setup structures to time functions
  struct timespec time1, time2, time3;
  double elapsed1 = 0.0, elapsed2 = 0.0;

  // set number of threads for parallel functions
  omp_set_num_threads(4);
  std::cout << std::setprecision(12);

  // get input parameters
  Parameters.density = atof(argv[1]);
  Parameters.Shells = atoi(argv[2]);
  Parameters.Pshells = atoi(argv[3]);
  Parameters.Nshells = atoi(argv[4]);
  Parameters.MBPT_Approx = atoi(argv[5]);
  Parameters.MBPT_Function = atoi(argv[6]);
  Build_Model_Space(Parameters, Space);
  
  clock_gettime(CLOCK_MONOTONIC, &time1); // time1 before channel function

  // setup channels for MBPT_Function == 2, 3, 4
  if(Parameters.MBPT_Function > 1){
    Setup_Channels_MBPT(Parameters, Space, Chan);
  }
  
  clock_gettime(CLOCK_MONOTONIC, &time2); // time2 after channel function, before MBPT function
  elapsed1 = (time2.tv_sec - time1.tv_sec);
  elapsed1 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;

  // Perform MBPT function according to MBPT_Approx and MBPT_Function
  if(Parameters.MBPT_Approx == 2){
    if(Parameters.MBPT_Function == 0){ MBPT2_0(Parameters, Space); }
    else if(Parameters.MBPT_Function == 1){ MBPT2_1(Parameters, Space); }
    else if(Parameters.MBPT_Function == 2){ MBPT2_2(Parameters, Space, Chan); }
    else if(Parameters.MBPT_Function == 3){ MBPT2_3(Parameters, Space, Chan); }
    else if(Parameters.MBPT_Function == 4){ MBPT2_4(Parameters, Space, Chan); }
  }
  else if(Parameters.MBPT_Approx == 3){
    if(Parameters.MBPT_Function == 0){ MBPT3_0(Parameters, Space); }
    else if(Parameters.MBPT_Function == 1){ MBPT3_1(Parameters, Space); }
    else if(Parameters.MBPT_Function == 2){ MBPT3_2(Parameters, Space, Chan); }
    else if(Parameters.MBPT_Function == 3){ MBPT3_3(Parameters, Space, Chan); }
    else if(Parameters.MBPT_Function == 4){ MBPT3_4(Parameters, Space, Chan); }
  }
  
  clock_gettime(CLOCK_MONOTONIC, &time3); // time3 after MBPT function
  elapsed2 = (time3.tv_sec - time2.tv_sec);
  elapsed2 += (time3.tv_nsec - time2.tv_nsec) / 1000000000.0;
  
  std::cout << "size = " << Space.indtot << ",  Approx = " << Parameters.MBPT_Approx << ",  Function = " << Parameters.MBPT_Function << ",  time1 = " << elapsed1 << ",  time2 = " << elapsed2 << std::endl << std::endl;
  
  return 0;
}
