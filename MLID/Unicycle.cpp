// Created by in July 2017, Codes inherited from Marcel Missura
#include "Unicycle.h"
#include "utility.h"
#include "fresnelnr.h"

using namespace std;

Unicycle::Unicycle() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
	v = 0.0;
	w = 0.0;
}

void Unicycle::clear() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
	v = 0.0;
	w = 0.0;
}

void Unicycle::set_state(double x, double y, double z, double v, double w) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->v = v;
	this->w = w;
}

void Unicycle::set_state(vec state) {
	this->x = state(0);
	this->y = state(1);
	this->z = state(2);
	this->v = state(3);
	this->w = state(4);
}

// Updates the car by time t using Euler integration with a timestep of 10^-6.
void Unicycle::simUpdateE6(double a, double b, double t)
{
	double timestep = 0.000001; // simulation time step

								//	double z0 = z;
								//	double v0 = v;
								//	double w0 = w;

	double time = 0;
	while (time <= t)
	{
		x += timestep * v * cos(z);
		y += timestep * v * sin(z);
		z += timestep * w;
		v += timestep * a;
		w += timestep * b;

		//		z = 0.5*b*time*time+w0*time+z0;
		//		v = a*time + v0;
		//		w = b*time + w0;

		time += timestep;
	}
	//z = normalize_angle(z);
}

// Updates the car analytically given the controls a, b and t.
void Unicycle::update(double a, double b, double t)
{
	double dx, dy, dz, dv, dw;

	double s = v;
	double c = w;
	double d = z;

	// The easy ones.
	dz = 0.5*b*t*t + c*t;
	dv = a*t;
	dw = b*t;

	// dx and dy

	// b = 0 and c = 0 case -> just linear acceleration
	if (fabs(b) < 0.001 && fabs(c) < 0.001)
	{
		double l = 0.5*a*t*t + s*t;
		dx = cos(d) * l;
		dy = sin(d) * l;
	}

	// b = 0 case -> formula singularity covered
	else if (fabs(b) < 0.001)
	{
		dx = a*(cos(c*t + d) - cos(d)) / (c*c) + ((a*t + s)*sin(c*t + d) - s*sin(d)) / c;
		dy = a*(sin(c*t + d) - sin(d)) / (c*c) - ((a*t + s)*cos(c*t + d) - s*cos(d)) / c;
	}

	else
	{
		// Use a mirroring technique if b is negative to avoid the negative sqrt().
		bool flipped = false;
		if (b < 0)
		{
			b = -b;
			c = -c;
			flipped = true;
		}

		double sb = sqrt(b);
		double pb15 = pow(b, 1.5);
		double gamma = cos(0.5*c*c / b - d);
		double sigma = sin(0.5*c*c / b - d);
		double c1, s1, c0, s0;
		frenel((c + b*t) / (sb*SPI), &s1, &c1);
		frenel(c / (sb*SPI), &s0, &c0);
		double C = c1 - c0;
		double S = s1 - s0;

		dx = SPI * (b*s - a*c) * (sigma*S + gamma*C) / pb15 + (a / b)*(sin(0.5*b*t*t + c*t + d) - sin(d));
		dy = SPI * (b*s - a*c) * (gamma*S - sigma*C) / pb15 - (a / b)*(cos(0.5*b*t*t + c*t + d) - cos(d));

		if (flipped)
		{
			double c2d = cos(2 * d);
			double s2d = sin(2 * d);
			double dxt = c2d*dx + s2d*dy;
			dy = s2d*dx - c2d*dy;
			dx = dxt;
		}
	}

	// Apply the state change and transform into world coordinates.
	x += dx;
	y += dy;
	z += dz;
	v += dv;
	w += dw;

	//w = normalize_angle(w);

	return;
}

void Unicycle::update(vec u) {
	update(u(0), u(1), u(2));
}

void Unicycle::update_3b(vec u_3b) {
	update(u_3b(0), u_3b(1), u_3b(2));
	update(u_3b(3), u_3b(4), u_3b(5));
	update(u_3b(6), u_3b(7), u_3b(8));
}

Mat<double> Unicycle::jacobian(double a, double b, double t)
{
	double dxda, dxdb, dxdt;
	double dyda, dydb, dydt;
	double dzda, dzdb, dzdt;
	double dvda, dvdb, dvdt;
	double dwda, dwdb, dwdt;

	double s = v;
	double c = w;

	// The easy ones.
	dxdt = (a*t + s) * cos(0.5*b*t*t + c*t);
	dydt = (a*t + s) * sin(0.5*b*t*t + c*t);
	dzda = 0.0;
	dzdb = 0.5*t*t;
	dzdt = b*t + c;
	dvda = t;
	dvdb = 0.0;
	dvdt = a;
	dwda = 0;
	dwdb = t;
	dwdt = b;

	// Difficult ones with zero conditions and b sign flip.

	// b = 0 and c = 0 case -> just linear acceleration
	if (fabs(b) < 0.00001 && fabs(c) < 0.00001)
	{
		dxda = 0.5*t*t;
		dyda = 0.0;

		// There is a singularity at b = 0.
		// Approximate the b gradient with an orthogonal unit vector.
		// This could be replaced by a more accurate "bridge".
		dxdb = 0.0;
		dydb = 1.0;
	}

	// b = 0 case
	else if (fabs(b) < 0.00001)
	{
		dxda = (c*t * sin(c*t) + cos(c*t) - 1) / (c*c);
		dyda = (sin(c*t) - c*t * cos(c*t)) / (c*c);

		// There is a singularity at b = 0.
		// Approximate the b gradient with an orthogonal unit vector.
		// This could be replaced by a more accurate "bridge".
		dxdb = -sin(c*t);
		dydb = cos(c*t);
	}

	else
	{
		double sgn = 1.0;
		if (b < 0)
		{
			b = -b;
			c = -c;
			sgn = -1.0;
		}

		double sb = sqrt(b);
		double pb15 = pow(b, 1.5);
		double pb25 = pow(b, 2.5);
		double pb35 = pow(b, 3.5);
		double gamma = cos(0.5*c*c / b);
		double sigma = sin(0.5*c*c / b);
		double kappa = cos(0.5*t*t*b + c*t);
		double zeta = sin(0.5*t*t*b + c*t);
		double c1, s1, c0, s0;
		frenel((c + b*t) / (sb*SPI), &s1, &c1);
		frenel(c / (sb*SPI), &s0, &c0);
		double C = c1 - c0;
		double S = s1 - s0;
		double dSdb = sin((b*t + c)*(b*t + c) / (2.0*b)) * (t / (SPI*sb) - (c + t*b) / (2.0*SPI*pb15)) + c*sigma / (2.0*SPI*pb15);
		double dCdb = cos((b*t + c)*(b*t + c) / (2.0*b)) * (t / (SPI*sb) - (c + t*b) / (2.0*SPI*pb15)) + c*gamma / (2.0*SPI*pb15);

		dxda = SPI*c / pb15 * (gamma*C - sigma*S) + zeta / b;
		dxdb = sgn * (SPI * (s / sb - a*c / pb15) * (sigma*dSdb + gamma*dCdb)
			+ SPI * (3.0*a*c / (2.0*pb25) - s / (2.0*pb15)) * (sigma*S + gamma*C)
			+ SPI * (a*c*c*c / (2.0*pb35) - s*c*c / (2.0*pb25)) * (gamma*S - sigma*C)
			+ a*t*t / (2 * b) * kappa - a / (b*b)*zeta);
		dyda = sgn * (SPI*c / pb15 * (sigma*C - gamma*S) - (kappa - 1.0) / b);
		dydb = SPI * (s / sb - a*c / pb15) * (gamma*dSdb - sigma*dCdb)
			+ SPI * (3.0*a*c / (2.0*pb25) - s / (2.0*pb15)) * (gamma*S - sigma*C)
			+ SPI * (-a*c*c*c / (2.0*pb35) + s*c*c / (2.0*pb25)) * (sigma*S + gamma*C)
			+ a*t*t / (2.0*b) * zeta + a / (b*b)*(kappa - 1.0);
	}


	// Assign values to Jacobian and rotate the x and y components by z to denormalize.
	Mat<double> J(5, 3);
	double cz = cos(z);
	double sz = sin(z);
	J(0, 0) = cz*dxda - sz*dyda;	J(0, 1) = cz*dxdb - sz*dydb;	J(0, 2) = cz*dxdt - sz*dydt;
	J(1, 0) = sz*dxda + cz*dyda;	J(1, 1) = sz*dxdb + cz*dydb;	J(1, 2) = sz*dxdt + cz*dydt;
	J(2, 0) = dzda;				J(2, 1) = dzdb;				J(2, 2) = dzdt;
	J(3, 0) = dvda;				J(3, 1) = dvdb;				J(3, 2) = dvdt;
	J(4, 0) = dwda;				J(4, 1) = dwdb;				J(4, 2) = dwdt;

	return J;
}

vec Unicycle::get_state_vec() {
	vec r(5);
	r(0) = x;
	r(1) = y;
	r(2) = z;
	r(3) = v;
	r(4) = w;
	return r;
}

// The second order Jacobian describes the effects of the first parameters on the second bang.
Mat<double> Unicycle::jacobian2(double a1, double b1, double t1, double a, double b, double t)
{
	double dxda1, dxdb1, dxdt1;
	double dyda1, dydb1, dydt1;
	double dzda1, dzdb1, dzdt1;
	double dvda1, dvdb1, dvdt1;
	double dwda1, dwdb1, dwdt1;

	double s = v;
	double c = w;
	double d = z;

	// b = 0 case. Singularity!
	if (fabs(b) < 0.00001)
	{
		double gamma = cos(c*t + d);
		double sigma = sin(c*t + d);
		double gamma0 = cos(d);
		double sigma0 = sin(d);

		if (fabs(c) < 0.00001)
		{
			dxda1 = cos(d)*t*t1;
			dyda1 = sin(d)*t*t1;

			// Approximate the b1 gradient with an orthogonal unit vector.
			// Replace this with a bridge.
			dxdb1 = -sin(c*t);
			dydb1 = cos(c*t);

			dxdt1 = cos(d)*a1*t;
			dydt1 = sin(d)*a1*t;
		}
		else
		{
			dxda1 = t1 / c * (sigma - sigma0);
			dyda1 = t1 / c * (gamma0 - gamma);

			// These are wrong.
			dxdb1 = t1*(t*(a*t + s)*gamma / c + sigma*(s - 2 * (a*t + s)) / (c*c) + 2 * a*(1 - gamma) / (c*c*c));
			dydb1 = t1*(t*(a*t + s)*sigma / c + (gamma - 1)*(2 * (a*t + s) - s) / (c*c) - 2 * a*sigma / (c*c*c));

			// These are wrong.
			dxdt1 = b1*(t*(a*t + s)*gamma / c + sigma*(s - 2 * (a*t + s)) / (c*c) + 2 * a*(1 - gamma) / (c*c*c)) + a1*sigma / c;
			dydt1 = b1*(t*(a*t + s)*sigma / c + (gamma - 1)*(2 * (a*t + s) - s) / (c*c) - 2 * a*sigma / (c*c*c)) + a1*(1 - gamma) / c;
		}

		dzda1 = 0.0;
		dvda1 = 0.0;
		dwda1 = 0.0;

		dzdb1 = t*t1;
		dvdb1 = 0.0;
		dwdb1 = 0.0;

		dzdt1 = t*b1;
		dvdt1 = 0.0;
		dwdt1 = 0.0;
	}

	else
	{
		// Use a mirroring technique if b is negative to avoid the negative sqrt().
		bool flipped = false;
		if (b < 0)
		{
			b = -b;
			b1 = -b1;
			c = -c;
			flipped = true;
		}

		double sb = sqrt(b);
		double pb15 = pow(b, 1.5);
		double sigma = sin(0.5*c*c / b - d);
		double gamma = cos(0.5*c*c / b - d);
		double zeta = sin(0.5*t*t*b + c*t + d);
		double kappa = cos(0.5*t*t*b + c*t + d);
		double s_ = sin((b*t + c)*(b*t + c) / (2.0*b));
		double c_ = cos((b*t + c)*(b*t + c) / (2.0*b));
		double s0_ = sin(c*c / (2.0*b));
		double c0_ = cos(c*c / (2.0*b));
		double sind = sin(d);
		double cosd = cos(d);
		double c1, s1, c0, s0;
		frenel((c + b*t) / (sb*SPI), &s1, &c1);
		frenel(c / (sb*SPI), &s0, &c0);
		double C = c1 - c0;
		double S = s1 - s0;
		double dSdb1 = (s_ - s0_) * t1 / (SPI*sb);
		double dCdb1 = (c_ - c0_) * t1 / (SPI*sb);
		double dSdt1 = (s_ - s0_) * b1 / (SPI*sb);
		double dCdt1 = (c_ - c0_) * b1 / (SPI*sb);

		dxda1 = t1 * SPI / sb * (sigma*S + gamma*C);
		dyda1 = t1 * SPI / sb * (gamma*S - sigma*C);
		dzda1 = 0.0;
		dvda1 = 0.0;
		dwda1 = 0.0;

		dxdb1 = SPI / pb15 * ((b*s - a*c) * (sigma*dSdb1 + gamma*dCdb1)
			+ (b*s - a*c) * (t1*c / b - 0.5*t1*t1) * (gamma*S - sigma*C)
			- t1 * a * (sigma*S + gamma*C))
			+ a / b * ((0.5*t1*t1 + t*t1)*kappa - 0.5*t1*t1*cosd);
		dydb1 = SPI / pb15 * ((b*s - a*c) * (gamma*dSdb1 - sigma*dCdb1)
			- (b*s - a*c) * (t1*c / b - 0.5*t1*t1) * (sigma*S + gamma*C)
			+ t1 * a * (sigma*C - gamma*S))
			+ a / b * ((0.5*t1*t1 + t*t1)*zeta - 0.5*t1*t1*sind);
		dzdb1 = t*t1;
		dvdb1 = 0.0;
		dwdb1 = 0.0;

		dxdt1 = SPI / pb15 * ((b*s - a*c) * (sigma*dSdt1 + gamma*dCdt1)
			+ (a1*b - a*b1) * (sigma*S + gamma*C)
			+ (b*s - a*c) * (b1*c / b - c) * (gamma*S - sigma*C))
			+ a / b * ((b1*t + c)*kappa - c*cosd);
		dydt1 = SPI / pb15 * ((b*s - a*c) * (gamma*dSdt1 - sigma*dCdt1)
			+ (a1*b - a*b1) * (gamma*S - sigma*C)
			- (b*s - a*c) * (b1*c / b - c) * (sigma*S + gamma*C))
			+ a / b * ((b1*t + c)*zeta - c*sind);
		dzdt1 = t*b1;
		dvdt1 = 0.0;
		dwdt1 = 0.0;

		if (flipped)
		{
			double dxt;
			double cos2d = cos(2 * d);
			double sin2d = sin(2 * d);

			dxt = cos2d*dxda1 + sin2d*dyda1;
			dyda1 = sin2d*dxda1 - cos2d*dyda1;
			dxda1 = dxt;

			dxt = cos2d*dxdb1 + sin2d*dydb1;
			dydb1 = -(sin2d*dxdb1 - cos2d*dydb1);
			dxdb1 = -dxt;

			dxt = cos2d*dxdt1 + sin2d*dydt1;
			dydt1 = (sin2d*dxdt1 - cos2d*dydt1);
			dxdt1 = dxt;
		}
	}


	// Assign values to the second order Jacobian.
	Mat<double> J(5, 3);

	J(0, 0) = dxda1; J(0, 1) = dxdb1; J(0, 2) = dxdt1;
	J(1, 0) = dyda1; J(1, 1) = dydb1; J(1, 2) = dydt1;
	J(2, 0) = dzda1; J(2, 1) = dzdb1; J(2, 2) = dzdt1;
	J(3, 0) = dvda1; J(3, 1) = dvdb1; J(3, 2) = dvdt1;
	J(4, 0) = dwda1; J(4, 1) = dwdb1; J(4, 2) = dwdt1;

	return J;
}

// The third order Jacobian describes the effects of the first parameters on the third link.
Mat<double> Unicycle::jacobian3(double a1, double b1, double t1, double a2, double b2, double t2, double a3, double b3, double t3)
{
	double dxda1, dxdb1, dxdt1;
	double dyda1, dydb1, dydt1;
	double dzda1, dzdb1, dzdt1;
	double dvda1, dvdb1, dvdt1;
	double dwda1, dwdb1, dwdt1;

	double s = v;
	double c = w;
	double d = z;
	double a = a3;
	double b = b3;
	double t = t3;

	// b = 0 case. Singularity!
	if (fabs(b) < 0.00001)
	{
		double gamma = cos(c*t + d);
		double sigma = sin(c*t + d);
		double gamma0 = cos(d);
		double sigma0 = sin(d);

		if (fabs(c) < 0.00001)
		{
			dxda1 = cos(d)*t*t1;
			dyda1 = sin(d)*t*t1;

			// Approximate the b1 gradient with an orthogonal unit vector.
			// Replace this with a bridge.
			dxdb1 = -sin(c*t);
			dydb1 = cos(c*t);

			dxdt1 = cos(d)*a1*t;
			dydt1 = sin(d)*a1*t;
		}
		else
		{
			dxda1 = t1 / c * (sigma - sigma0);
			dyda1 = t1 / c * (gamma0 - gamma);

			// These are wrong.
			dxdb1 = t1*(t*(a*t + s)*gamma / c + sigma*(s - 2 * (a*t + s)) / (c*c) + 2 * a*(1 - gamma) / (c*c*c));
			dydb1 = t1*(t*(a*t + s)*sigma / c + (gamma - 1)*(2 * (a*t + s) - s) / (c*c) - 2 * a*sigma / (c*c*c));

			// These are wrong.
			dxdt1 = b1*(t*(a*t + s)*gamma / c + sigma*(s - 2 * (a*t + s)) / (c*c) + 2 * a*(1 - gamma) / (c*c*c)) + a1*sigma / c;
			dydt1 = b1*(t*(a*t + s)*sigma / c + (gamma - 1)*(2 * (a*t + s) - s) / (c*c) - 2 * a*sigma / (c*c*c)) + a1*(1 - gamma) / c;
		}

		dzda1 = 0.0;
		dvda1 = 0.0;
		dwda1 = 0.0;

		dzdb1 = t*t1;
		dvdb1 = 0.0;
		dwdb1 = 0.0;

		dzdt1 = t*b1;
		dvdt1 = 0.0;
		dwdt1 = 0.0;
	}

	else
	{
		// Use a mirroring technique if b is negative to avoid the negative sqrt().
		bool flipped = false;
		if (b < 0)
		{
			b = -b;
			b1 = -b1;
			b2 = -b2;
			c = -c;
			flipped = true;
		}

		double sb = sqrt(b);
		double pb15 = pow(b, 1.5);
		double sigma = sin(0.5*c*c / b - d);
		double gamma = cos(0.5*c*c / b - d);
		double zeta = sin(0.5*t*t*b + c*t + d);
		double kappa = cos(0.5*t*t*b + c*t + d);
		double sind = sin(d);
		double cosd = cos(d);
		double c1, s1, c0, s0;
		frenel((c + b*t) / (sb*SPI), &s1, &c1);
		frenel(c / (sb*SPI), &s0, &c0);
		double C = c1 - c0;
		double S = s1 - s0;
		double s_ = sin((b*t + c)*(b*t + c) / (2.0*b));
		double c_ = cos((b*t + c)*(b*t + c) / (2.0*b));
		double s0_ = sin(c*c / (2.0*b));
		double c0_ = cos(c*c / (2.0*b));
		double dSdb1 = (s_ - s0_) * t1 / (SPI*sb);
		double dCdb1 = (c_ - c0_) * t1 / (SPI*sb);
		double dSdt1 = (s_ - s0_) * b1 / (SPI*sb);
		double dCdt1 = (c_ - c0_) * b1 / (SPI*sb);

		dxda1 = t1 * SPI / sb * (sigma*S + gamma*C);
		dyda1 = t1 * SPI / sb * (gamma*S - sigma*C);
		dzda1 = 0.0;
		dvda1 = 0.0;
		dwda1 = 0.0;

		dxdb1 = SPI / pb15 * ((b*s - a*c) * (sigma*dSdb1 + gamma*dCdb1)
			+ (b*s - a*c) * (t1*c / b - t1*t2 - 0.5*t1*t1) * (gamma*S - sigma*C)
			- t1 * a * (sigma*S + gamma*C))
			+ a / b * ((0.5*t1*t1 + t1*t2 + t*t1)*kappa - (0.5*t1*t1 + t1*t2)*cosd);
		dydb1 = SPI / pb15 * ((b*s - a*c) * (gamma*dSdb1 - sigma*dCdb1)
			- (b*s - a*c) * (t1*c / b - t1*t2 - 0.5*t1*t1) * (sigma*S + gamma*C)
			+ t1 * a * (sigma*C - gamma*S))
			+ a / b * ((0.5*t1*t1 + t1*t2 + t*t1)*zeta - (0.5*t1*t1 + t1*t2)*sind);
		dzdb1 = t*t1;
		dvdb1 = 0.0;
		dwdb1 = 0.0;

		dxdt1 = SPI / pb15 * ((b*s - a*c) * (sigma*dSdt1 + gamma*dCdt1)
			+ (a1*b - a*b1) * (sigma*S + gamma*C)
			+ (b*s - a*c) * (b1*c / b - b1*t2 - (c - b2*t2)) * (gamma*S - sigma*C))
			+ a / b * ((b1*t2 + b1*t + (c - b2*t2))*kappa - (b1*t2 + (c - b2*t2))*cosd);
		dydt1 = SPI / pb15 * ((b*s - a*c) * (gamma*dSdt1 - sigma*dCdt1)
			+ (a1*b - a*b1) * (gamma*S - sigma*C)
			- (b*s - a*c) * (b1*c / b - b1*t2 - (c - b2*t2)) * (sigma*S + gamma*C))
			+ a / b * ((b1*t2 + b1*t + (c - b2*t2))*zeta - (b1*t2 + (c - b2*t2))*sind);
		dzdt1 = t*b1;
		dvdt1 = 0.0;
		dwdt1 = 0.0;

		if (flipped)
		{
			double dxt;
			double cos2d = cos(2 * d);
			double sin2d = sin(2 * d);

			dxt = cos2d*dxda1 + sin2d*dyda1;
			dyda1 = sin2d*dxda1 - cos2d*dyda1;
			dxda1 = dxt;

			dxt = cos2d*dxdb1 + sin2d*dydb1;
			dydb1 = -(sin2d*dxdb1 - cos2d*dydb1);
			dxdb1 = -dxt;

			dxt = cos2d*dxdt1 + sin2d*dydt1;
			dydt1 = (sin2d*dxdt1 - cos2d*dydt1);
			dxdt1 = dxt;
		}
	}

	// Assign values to the third order Jacobian.
	Mat<double> J(5, 3);

	J(0, 0) = dxda1; J(0, 1) = dxdb1; J(0, 2) = dxdt1;
	J(1, 0) = dyda1; J(1, 1) = dydb1; J(1, 2) = dydt1;
	J(2, 0) = dzda1; J(2, 1) = dzdb1; J(2, 2) = dzdt1;
	J(3, 0) = dvda1; J(3, 1) = dvdb1; J(3, 2) = dvdt1;
	J(4, 0) = dwda1; J(4, 1) = dwdb1; J(4, 2) = dwdt1;

	return J;
}

// The complete Jacobian matrix of a three link trajectory.
Mat<double> Unicycle::jacobian(double a1, double b1, double t1, double a2, double b2, double t2, double a3, double b3, double t3)
{
	Mat<double> J(5, 3);
	double dxda1, dxdb1, dxdt1, dxda2, dxdb2, dxdt2, dxda3, dxdb3, dxdt3;
	double dyda1, dydb1, dydt1, dyda2, dydb2, dydt2, dyda3, dydb3, dydt3;
	double dzda1, dzdb1, dzdt1, dzda2, dzdb2, dzdt2, dzda3, dzdb3, dzdt3;
	double dvda1, dvdb1, dvdt1, dvda2, dvdb2, dvdt2, dvda3, dvdb3, dvdt3;
	double dwda1, dwdb1, dwdt1, dwda2, dwdb2, dwdt2, dwda3, dwdb3, dwdt3;

	Unicycle car(x, y, z, v, w);

	// The first half of the Jacobian is the one link Jacobian of the first state
	// plus the "second order" Jacobian of the second state
	// plus the "third order" Jacobian of the third state.
	J = jacobian(a1, b1, t1);

	car.z += car.w*t1 + 0.5*b1*t1*t1;
	car.v += a1*t1;
	car.w += b1*t1;
	J += car.jacobian2(a1, b1, t1, a2, b2, t2);

	car.z += car.w*t2 + 0.5*b2*t2*t2;
	car.v += a2*t2;
	car.w += b2*t2;
	J += car.jacobian3(a1, b1, t1, a2, b2, t2, a3, b3, t3);

	dxda1 = J(0, 0);	dxdb1 = J(0, 1);	dxdt1 = J(0, 2);
	dyda1 = J(1, 0);	dydb1 = J(1, 1);	dydt1 = J(1, 2);
	dzda1 = J(2, 0);	dzdb1 = J(2, 1);	dzdt1 = J(2, 2);
	dvda1 = J(3, 0);	dvdb1 = J(3, 1);	dvdt1 = J(3, 2);
	dwda1 = J(4, 0);	dwdb1 = J(4, 1);	dwdt1 = J(4, 2);


	// The second half of the Jacobian is the one link Jacobian of the second state
	// plus the second order Jacobian from the third state.
	car = Unicycle(x, y, z, v, w);
	car.z += car.w*t1 + 0.5*b1*t1*t1;
	car.v += a1*t1;
	car.w += b1*t1;
	J = car.jacobian(a2, b2, t2);

	car.z += car.w*t2 + 0.5*b2*t2*t2;
	car.v += a2*t2;
	car.w += b2*t2;
	J += car.jacobian2(a2, b2, t2, a3, b3, t3);

	dxda2 = J(0, 0);	dxdb2 = J(0, 1);	dxdt2 = J(0, 2);
	dyda2 = J(1, 0);	dydb2 = J(1, 1);	dydt2 = J(1, 2);
	dzda2 = J(2, 0);	dzdb2 = J(2, 1);	dzdt2 = J(2, 2);
	dvda2 = J(3, 0);	dvdb2 = J(3, 1);	dvdt2 = J(3, 2);
	dwda2 = J(4, 0);	dwdb2 = J(4, 1);	dwdt2 = J(4, 2);


	// The third half of the three link Jacobian is the one link Jacobian of the third state.
	J = car.jacobian(a3, b3, t3);
	dxda3 = J(0, 0);	dxdb3 = J(0, 1);	dxdt3 = J(0, 2);
	dyda3 = J(1, 0);	dydb3 = J(1, 1);	dydt3 = J(1, 2);
	dzda3 = J(2, 0);	dzdb3 = J(2, 1);	dzdt3 = J(2, 2);
	dvda3 = J(3, 0);	dvdb3 = J(3, 1);	dvdt3 = J(3, 2);
	dwda3 = J(4, 0);	dwdb3 = J(4, 1);	dwdt3 = J(4, 2);


	// Assign values to the three link Jacobian.
	Mat<double> J3(5, 9);
	J3(0, 0) = dxda1; J3(0, 1) = dxdb1; J3(0, 2) = dxdt1; J3(0, 3) = dxda2; J3(0, 4) = dxdb2; J3(0, 5) = dxdt2; J3(0, 6) = dxda3; J3(0, 7) = dxdb3; J3(0, 8) = dxdt3;
	J3(1, 0) = dyda1; J3(1, 1) = dydb1; J3(1, 2) = dydt1; J3(1, 3) = dyda2; J3(1, 4) = dydb2; J3(1, 5) = dydt2; J3(1, 6) = dyda3; J3(1, 7) = dydb3; J3(1, 8) = dydt3;
	J3(2, 0) = dzda1; J3(2, 1) = dzdb1; J3(2, 2) = dzdt1; J3(2, 3) = dzda2; J3(2, 4) = dzdb2; J3(2, 5) = dzdt2; J3(2, 6) = dzda3; J3(2, 7) = dzdb3; J3(2, 8) = dzdt3;
	J3(3, 0) = dvda1; J3(3, 1) = dvdb1; J3(3, 2) = dvdt1; J3(3, 3) = dvda2; J3(3, 4) = dvdb2; J3(3, 5) = dvdt2; J3(3, 6) = dvda3; J3(3, 7) = dvdb3; J3(3, 8) = dvdt3;
	J3(4, 0) = dwda1; J3(4, 1) = dwdb1; J3(4, 2) = dwdt1; J3(4, 3) = dwda2; J3(4, 4) = dwdb2; J3(4, 5) = dwdt2; J3(4, 6) = dwda3; J3(4, 7) = dwdb3; J3(4, 8) = dwdt3;

	return J3;
}