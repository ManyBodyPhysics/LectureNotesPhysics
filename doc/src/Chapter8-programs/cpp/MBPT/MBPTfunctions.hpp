#ifndef CCFUNCTIONS_H
#define CCFUNCTIONS_H

#include "parameters.hpp"
#include "channels.hpp"
#include "modelspace.hpp"

void MBPT3_0(const Input_Parameters &Parameters, const Model_Space &Space);
void MBPT3_1(const Input_Parameters &Parameters, const Model_Space &Space);
void MBPT3_2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT3_3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT3_4(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);

#endif
