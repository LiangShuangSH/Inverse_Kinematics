#include "utility.h"

double normalize_angle(double angle) {
	double new_angle = fmod(angle + PI, 2*PI);
	if (new_angle < 0) {
		new_angle += 2 * PI;
	}
	return new_angle - PI;
}

double euclidean_distance(vec a, vec b) {
	double result = 0.0;
	if (a.size() != b.size()) {
		cout << "Warning::Euclidean_Distance()::SIZE_NOT_EQUAL" << endl;
		result = 0.0;
	}
	else {
		for (int i = 0; i < a.size(); i++) {
			result += pow(a(i) - b(i), 2);
		}
		result = sqrt(result);
	}
	return result;
}

vec normalize_vector(vec v, int dim, double scale) {
	double IVI = 0.0;
	for (int i = 0; i < dim; i++) {
		IVI += pow(v(i), 2);
	}
	IVI = sqrt(IVI);
	vec result = scale * v / IVI;

	return result;
}