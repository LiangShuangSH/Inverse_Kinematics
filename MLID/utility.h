#include <math.h>
#include "globals.h"
#include <armadillo>

using namespace arma;

double normalize_angle(double angle);
double euclidean_distance(vec a, vec b);
vec normalize_vector(vec v, int dim, double scale);