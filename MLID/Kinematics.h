#ifndef Kinematics_H
#define Kinematics_H

#include <iostream>
#include "Unicycle.h"
#include "globals.h"

const double A_MAX = 3.0;
const double A_MIN = -3.0;
const double B_MAX = PI;
const double B_MIN = -PI;
const double T_MAX = 10.0;
const double T_MIN = 0.0;
const double V_MAX = 120000.0/3600.0;
const double W_MAX = PI;

vec FK(vec s0, vec u);
vec FK_3b(vec s0, vec u);
vec IK(vec sq, vec ui, vec s0);
vec IK_3b(vec sq, vec uk, vec s0);
vec random_uk(int dim);
vec random_s();

#endif