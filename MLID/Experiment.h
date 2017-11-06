#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include <iostream>
#include <fstream>
#include "Unicycle.h"
#include "utility.h"
#include "Kinematics.h"
#include "Sampler.h"

using namespace std;

void Experiment();
void Experiment2();
bool check_conditions(vec s0, vec uk);
void sample(int num);

#endif /* EXPERIMENT_H_ */