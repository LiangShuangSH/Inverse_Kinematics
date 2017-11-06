#ifndef Sampler_H
#define Sampler_H

#include "Kinematics.h"
#include <fstream>

using namespace std;

void sampling(vec s0, double itv);
void uniform_sample_uspace(vec s0, Mat<double> &uks, Mat<double> &sks, double itv);
void write_data(Mat<double> uks, Mat<double> sks, bool overwrite);

#endif