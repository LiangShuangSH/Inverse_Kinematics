// Created by in July 2017, Codes inherited from Marcel Missura

#ifndef Unicycle_H
#define Unicycle_H

#include <math.h>
#include <armadillo>

using namespace arma;

class Unicycle {
public:
	double x;
	double y;
	double z;
	double v;
	double w;

	Unicycle();
	Unicycle(double xx, double yy, double zz, double vv, double ww) { x = xx; y = yy; z = zz; v = vv; w = ww; };
	~Unicycle() {};

	void clear();
	void set_state(double x, double y, double z, double v, double w);
	void set_state(vec state);
	void simUpdateE6(double a, double b, double t);
	void update(double a, double b, double t);
	void update(vec u);
	void update_3b(vec u_3b);
	Mat<double> jacobian(double a, double b, double t);
	Mat<double> jacobian2(double a1, double b1, double t1, double a, double b, double t);
	Mat<double> jacobian3(double a1, double b1, double t1, double a2, double b2, double t2, double a3, double b3, double t3);
	Mat<double> jacobian(double a1, double b1, double t1, double a2, double b2, double t2, double a3, double b3, double t3);
	vec get_state_vec();
};

#endif