// Created by in July 2017, Codes inherited from Marcel Missura
#ifndef GLOBS_H_
#define GLOBS_H_
#include <math.h>

const double G = 9.81;
const double PI = 3.1415926535897932384626433832795;
const double PII = 2.0*PI;
const double PI2 = 1.5707963267948965579989817342721;
const double PI4 = 0.78539816339744830961566084581988;
const double SPI = 1.7724538509055160272981674833411; // sqrt of pi
const double EPSILON = 1.0E-6;
const int MAX_POLYGON_VERTEXES = 16;

// Returns the sign of the argument, +1 if the argument is positive,
// -1 if the argument is negative, and 0 if the argument is zero.
template <typename T>
inline T sgn0(const T a) { return (a == 0 ? 0 : a < 0 ? -1 : 1); }

// Returns the sign of the argument, +1 if the argument is positive,
// and -1 if the argument is negative.
template <typename T>
inline T sgn(const T a) { return (a < 0 ? -1 : 1); }

// Maps the argument to an angle expressed in radians between -PI and PI.
inline double picut(double x) { return  x < 0 ? fmod(x - PI, PII) + PI : fmod(x + PI, PII) - PI; }

// Maps the argument to an angle expressed in radians between 0 and 2*PI.
inline double pi2cut(double x) { return fmod(fabs(x), PII); }

// A better version of modulo (%) that wraps rather than retuning a number < 0.
inline int mod(int x, int y) { return x >= 0 ? x%y : (x%y) + y; }

#endif /* GLOBS_H_ */
