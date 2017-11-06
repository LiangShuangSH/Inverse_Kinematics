#include <iostream>
#include <fstream>
#include "Unicycle.h"
#include <time.h>
#include <stdlib.h>
#include "globals.h"
#include <vector>
#include "utility.h"
#include "Kinematics.h"
#include "Sampler.h"
#include "Experiment.h"

using namespace std;

Mat<double> get_error_grad(vec s0, vec uk);
Mat<double> get_error_grad_3b(vec s0, vec uk);
Mat<double> rand_kernel_test(int unum);
Mat<double> kernel_test(int unum, double a, double b, double t);
Mat<double> rand_kernel_test_3b(int unum);
Mat<double> kernel_test_3b(int unum, vec uk);
vec get_state_error(vec eu, vec u);
Mat<double> rand_generator(int num, vec s0, double amin, double amax, double bmin, double bmax, double tmax);
void rand_kernel_ik(vec sq, vec s0);
vec gd_uk(vec s0, vec sq, vec uk0, double lr);
vec deu_duk(vec s0, vec sq, vec uk);
Mat<double> des_du(vec s0, vec uk);
double get_loss(vec s0, vec sq, vec uk);
void kernel_state_test(vec s0, vec uk, bool thrb);

int main() {
	
	srand(time(NULL));
	//srand(47);

	//Experiment();

	vec s0(5);
	s0.zeros();

	vec ut(9);
	ut << 1.02978 << -0.0143787 << 2.07238 << 1.90052 << 0.0882537 << 0.864398 << 2.17105 << -1.6757 << 0.24804 << endr;
	
	vec u2(9);
	u2 << 1.18854 << -0.00709698 << 1.97316 << 1.89831 << -0.0451131 << 0.820281 << 2.18084 << -1.65262 << 0.189247 << endr;

	vec u_error = u2 - ut;

	cout << "u_error: " << endl;
	u_error.print();
	cout << "" << endl;

	Unicycle uc;
	uc.set_state(s0);

	Mat<double> j = uc.jacobian(1.0376, 0.00508147, 1.98614, 1.93246, 0.10882, 0.85757, 2.18772, -1.66413, 0.233161);
	Mat<double> ji = pinv(j);
	cout << "J^-1: " << endl;
	ji.print();

	cout << "" << endl;

	Mat<double> e = j * u_error;
	cout << "e: " << endl;
	e.print();


	// Test Jacobian for 3 links
	/*
	Unicycle uc;
	uc.set_state(s0);
	vec uk1(9);
	uk1.fill(1.0);

	vec uk2 = uk1;
	uk2(1) += 0.1;
	vec s1 = FK_3b(s0, uk2);
	vec s2 = FK_3b(s0, uk1) 
		+ uc.jacobian(uk1(0), uk1(1), uk1(2),
		uk1(3), uk1(4), uk1(5),
		uk1(6), uk1(7), uk1(8)) * (uk2 - uk1);

	cout << "v: " << s0(3) << ", w: " << s0(4) << endl;

	cout << "s1: " << endl;
	s1.print();
	cout << "" << endl;
	cout << "s2: " << endl;
	s2.print();
	cout << "" << endl;

	vec s_error = s2 - s1;
	cout << "s_error: " << endl;
	s_error.print();
	*/

	// Inverse Kinematics Test
	/*
	vec uk1(9);

	while (true) {
		uk1 = random_uk(9);

		double angle1 = 0.0;
		double angle2 = 0.0;
		double angle3 = 0.0;

		double b1 = uk1(1);
		double b2 = uk1(4);
		double b3 = uk1(7);
		double t1 = uk1(2);
		double t2 = uk1(5);
		double t3 = uk1(8);

		double w0 = s0(4);
		double w1 = w0 + b1 * t1;
		double w2 = w1 + b2 * t2;

		angle1 = w0 * t1 + 0.5 * b1 * t1 * t1;
		angle2 = w1 * t2 + 0.5 * b2 * t2 * t2;
		angle3 = w2 * t3 + 0.5 * b3 * t3 * t3;

		if (
			abs(angle1) < PII
			&& abs(angle2) < PII
			&& abs(angle3) < PII
			&& abs(angle1 + angle2) < PII
			&& abs(angle2 + angle3) < PII
			&& abs(angle1 + angle2 + angle3) < PII
			) {
			break;
		}
	}

	vec sk1 = FK_3b(s0, uk1);

	vec uk2 = uk1;
	int idx = 2;
	double delta = 0.05;
	uk2(idx) += delta;

	vec sk2 = FK_3b(s0, uk2);

	vec uk2_i = IK_3b(sk2, uk1, s0);
	vec sk2_i = FK_3b(s0, uk2_i);

	double error = euclidean_distance(sk2_i, sk2);

	cout << "uk1 is " << endl;
	uk1.print();
	cout << "" << endl;

	cout << "uk(" << idx << ")" << " += " << delta << endl;
	cout << "" << endl;

	//cout << "¦¤x = " << sk2(0) - sk1(0) << endl;
	//cout << "¦¤y = " << sk2(1) - sk1(1) << endl;
	//cout << "¦¤theta = " << sk2(2) - sk1(2) << endl;
	//cout << "¦¤v = " << sk2(3) - sk1(3) << endl;
	//cout << "¦¤w = " << sk2(4) - sk1(4) << endl;
	//cout << "" << endl;

	cout << "uk2 inverse: " << endl;
	uk2_i.print();
	cout << "" << endl;

	cout << "error = " << error << endl;
	*/

	return 0;
}

// Call srand first!
Mat<double> rand_kernel_test(int unum) {

	Mat<double> record(3 + 3, 1 + unum);
	record.zeros();

	// Generate an arbitary kernel (s0, u0)
	vec s0(5);
	s0.zeros();
	vec u0(3);
	u0(0) = A_MIN + ((double)rand() / RAND_MAX) * (A_MAX - A_MIN);
	u0(1) = B_MIN + ((double)rand() / RAND_MAX) * (B_MAX - B_MIN);
	u0(2) = 0.0 + ((double)rand() / RAND_MAX) * (T_MAX - 0.0);
	for (int i = 0; i < 3; i++) {
		record(i, 0) = u0(i);
	}

	// Calculate the error gradient of the kernel
	Mat<double> error_grad(3, 3);
	error_grad = get_error_grad(s0, u0);

	// Generate random u around u0, record the error
	double a_rl = 0.3;
	double b_rl = 0.3;
	double t_rl = 0.3;
	for (int i = 1; i < unum+1; i++) {
		vec u(3);
		u(0) = u0(0) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(1) = u0(1) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(2) = u0(2) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		if (u(2) < 0.0) { u(2) = 0.0; }
		// Record u
		for (int j = 0; j < 3; j++) {
			record(j, i) = u(j);
		}
		// Record error of u
		vec error = error_grad * (u - u0);
		for (int j = 0; j < 3; j++) {
			record(j + 3, i) = error(j);
		}
	}

	return record;
}

Mat<double> rand_kernel_test_3b(int unum) {
	Mat<double> record(9 + 9, 1 + unum);
	record.zeros();

	// Generate an arbitary kernel (s0, u0)
	vec s0(5);
	s0.zeros();
	vec u0(9);
	u0(0) = A_MIN + ((double)rand() / RAND_MAX) * (A_MAX - A_MIN);
	u0(1) = B_MIN + ((double)rand() / RAND_MAX) * (B_MAX - B_MIN);
	u0(2) = 0.0 + ((double)rand() / RAND_MAX) * (T_MAX - 0.0);
	u0(3) = A_MIN + ((double)rand() / RAND_MAX) * (A_MAX - A_MIN);
	u0(4) = B_MIN + ((double)rand() / RAND_MAX) * (B_MAX - B_MIN);
	u0(5) = 0.0 + ((double)rand() / RAND_MAX) * (T_MAX - 0.0);
	u0(6) = A_MIN + ((double)rand() / RAND_MAX) * (A_MAX - A_MIN);
	u0(7) = B_MIN + ((double)rand() / RAND_MAX) * (B_MAX - B_MIN);
	u0(8) = 0.0 + ((double)rand() / RAND_MAX) * (T_MAX - 0.0);
	for (int i = 0; i < 9; i++) {
		record(i, 0) = u0(i);
	}

	// Calculate the error gradient of the kernel
	Mat<double> error_grad(9, 9);
	error_grad = get_error_grad_3b(s0, u0);

	// Generate random u around u0, record the error
	double a_rl = 0.3;
	double b_rl = 0.3;
	double t_rl = 0.3;
	for (int i = 1; i < unum + 1; i++) {
		vec u(9);
		u(0) = u0(0) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(1) = u0(1) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(2) = u0(2) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		u(3) = u0(3) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(4) = u0(4) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(5) = u0(5) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		u(6) = u0(6) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(7) = u0(7) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(8) = u0(8) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		if (u(2) < 0.0) { u(2) = 0.0; }
		if (u(5) < 0.0) { u(5) = 0.0; }
		if (u(8) < 0.0) { u(8) = 0.0; }
		// Record u
		for (int j = 0; j < 9; j++) {
			record(j, i) = u(j);
		}
		// Record error of u
		vec error = error_grad * (u - u0);
		for (int j = 0; j < 9; j++) {
			record(j + 9, i) = error(j);
		}
	}
	return record;
}

Mat<double> kernel_test(int unum, double a, double b, double t) {
	Mat<double> record(3 + 3, 1 + unum);
	record.zeros();

	// Generate an arbitary kernel (s0, u0)
	vec s0(5);
	s0.zeros();
	vec u0(3);
	u0(0) = a;
	u0(1) = b;
	u0(2) = t;
	for (int i = 0; i < 3; i++) {
		record(i, 0) = u0(i);
	}

	// Calculate the error gradient of the kernel
	Mat<double> error_grad(3, 3);
	error_grad = get_error_grad(s0, u0);

	// Generate random u around u0, record the error
	double a_rl = 10.0;
	double b_rl = 10.0;
	double t_rl = 0.1;
	for (int i = 1; i < unum + 1; i++) {
		vec u(3);
		u(0) = u0(0) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(1) = u0(1) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(2) = u0(2) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		if (u(2) < 0.0) { u(2) = 0.0; }
		// Record u
		for (int j = 0; j < 3; j++) {
			record(j, i) = u(j);
		}
		// Record error of u
		vec error = error_grad * (u - u0);
		for (int j = 0; j < 3; j++) {
			record(j + 3, i) = error(j);
		}
	}

	return record;
}

Mat<double> kernel_test_3b(int unum, vec uk) {
	Mat<double> record(9 + 9, 1 + unum);
	record.zeros();

	// Generate an arbitary kernel (s0, u0)
	vec s0(5);
	s0.zeros();
	vec u0(9);
	u0 = uk;
	for (int i = 0; i < 9; i++) {
		record(i, 0) = u0(i);
	}

	// Calculate the error gradient of the kernel
	Mat<double> error_grad(9, 9);
	error_grad = get_error_grad_3b(s0, u0);

	// Generate random u around u0, record the error
	double a_rl = 1.0;
	double b_rl = 1.0;
	double t_rl = 1.0;
	for (int i = 1; i < unum + 1; i++) {
		vec u(9);
		u(0) = u0(0) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(1) = u0(1) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(2) = u0(2) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		u(3) = u0(3) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(4) = u0(4) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(5) = u0(5) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		u(6) = u0(6) + (-a_rl) + ((double)rand() / RAND_MAX) * (2 * a_rl);
		u(7) = u0(7) + (-b_rl) + ((double)rand() / RAND_MAX) * (2 * b_rl);
		u(8) = u0(8) + (-t_rl) + ((double)rand() / RAND_MAX) * (2 * t_rl);
		if (u(2) < 0.0) { u(2) = 0.0; }
		if (u(5) < 0.0) { u(5) = 0.0; }
		if (u(8) < 0.0) { u(8) = 0.0; }
		// Record u
		for (int j = 0; j < 9; j++) {
			record(j, i) = u(j);
		}
		// Record error of u
		vec error = error_grad * (u - u0);
		for (int j = 0; j < 9; j++) {
			record(j + 9, i) = error(j);
		}
	}
	return record;
}

Mat<double> get_error_grad(vec s0, vec uk) {
	Mat<double> error_grad(3, 3);
	double small_value = 0.000001;

	Unicycle uc;
	vec u2(3);
	vec u3(3);
	for (int i = 0; i < 3; i++) {
		u2(0) = uk(0); u2(1) = uk(1); u2(2) = uk(2);
		u2(i) += small_value;

		uc.set_state(s0);
		uc.update(u2(0), u2(1), u2(2));

		u3 = IK(uc.get_state_vec(), uk, s0);

		vec partial_deri = (u3 - u2) / small_value;
		error_grad.col(i) = partial_deri;
	}

	return error_grad;
}

Mat<double> get_error_grad_3b(vec s0, vec uk) {
	Mat<double> error_grad(9, 9);
	double small_value = 0.01;

	Unicycle uc;
	vec u2(9);
	vec u3(9);

	for (int i = 0; i < 9; i++) {
		u2 = uk;
		u2(i) += small_value;

		uc.set_state(s0);
		uc.update_3b(u2);

		u3 = IK_3b(uc.get_state_vec(), uk, s0);

		vec partial_deri = (u3 - u2) / small_value;
		error_grad.col(i) = partial_deri;
	}
	return error_grad;
}

// eu: error of u
// u:  the u returned by inverse function
vec get_state_error(vec eu, vec u) {
	vec es(5);

	Unicycle uc1;
	Unicycle uc2;

	uc1.update(eu + u);
	uc2.update(u);

	es = uc1.get_state_vec() - uc2.get_state_vec();

	return es;
}

void rand_kernel_ik(vec sq, vec s0) {
	int k_num = 100000;

	vector<vec> u_nice;
	vector<vec> s_nice;

	// Generate random data points
	Mat<double> data = rand_generator(k_num, s0, A_MIN, A_MAX, B_MIN, B_MAX, T_MAX);

	for (int i = 0; i < k_num; i++) {
		vec uq = IK_3b(sq, data.col(i), s0);
		Unicycle uc;
		uc.set_state(s0);
		uc.update_3b(uq);
		vec sr = uc.get_state_vec();
		vec sd = sr - sq;
		if (abs(sd(0)) < 0.1
			&& abs(sd(1)) < 0.1
			&& abs(sd(2)) < 0.3
			&& abs(sd(3)) < 1.0
			&& abs(sd(4)) < 1.0) {
			u_nice.push_back(data.col(i));
			s_nice.push_back(sr);
		}
	}

	for (vec s:s_nice) {
		s.print();
		cout << "" << endl;
	}
}

Mat<double> rand_generator(int num, vec s0, double amin, double amax, double bmin, double bmax, double tmax) {
	Mat<double> data(9, num);

	for (int i = 0; i < num; i++) {
		vec u(9);
		u(0) = amin + ((double)rand() / RAND_MAX) * (amax - amin);
		u(1) = bmin + ((double)rand() / RAND_MAX) * (bmax - bmin);
		u(2) = 0.0 + ((double)rand() / RAND_MAX) * (tmax - 0.0);
		u(3) = amin + ((double)rand() / RAND_MAX) * (amax - amin);
		u(4) = bmin + ((double)rand() / RAND_MAX) * (bmax - bmin);
		u(5) = 0.0 + ((double)rand() / RAND_MAX) * (tmax - 0.0);
		u(6) = amin + ((double)rand() / RAND_MAX) * (amax - amin);
		u(7) = bmin + ((double)rand() / RAND_MAX) * (bmax - bmin);
		u(8) = 0.0 + ((double)rand() / RAND_MAX) * (tmax - 0.0);

		// Gimme a nice track
		double delta_z = 0.0;
		// delta_z = w0*t1 + 0.5*b1*t1^2 + w1*t2 + 0.5*b2*t2^2 + w2*t3 + 0.5*b3*t3^2
		double b1 = u(1); double b2 = u(4); double b3 = u(7);
		double t1 = u(2); double t2 = u(5); double t3 = u(8);
		double w0 = s0(4); double w1 = w0 + b1*t1; double w2 = w1 + b2*t2;
		delta_z = w0*t1 + 0.5*b1*pow(t1, 2) + w1*t2 + 0.5*b2*pow(t2, 2) + w2*t3 + 0.5*b3*pow(t3, 2);
		if (abs(delta_z) > PII) {
			i--;
			continue;
		}
		else {
			data.col(i) = u;
		}
	}
	return data;
}

vec gd_uk(vec s0, vec sq, vec uk0, double lr) {
	vec uk = uk0;
	vec uk2 = uk;
	bool converge = false;
	vec grad(9);

	int i = 0;
	while (!converge) {
		grad = deu_duk(s0, sq, uk);
		uk2 = uk - grad * lr;
		// Check converge
		/*
		vec ud = uk2 - uk;
		double diff = 0.0;
		for (int i = 0; i < 9; i++) {
			diff += pow(ud(i), 2);
		}
		if (diff < 0.0000000001) {
			converge = true;
		}
		*/
		
		
		i++;
		if (i == 30000) {
			converge = true;
		}
		
		
		/*
		if (get_loss(s0, sq, uk) < 0.02) {
			converge = true;
		}
		*/

		uk = uk2;
		uk.print();
		cout << "" << endl;
	}
	return uk;
}

vec deu_duk(vec s0, vec sq, vec uk) {
	vec deuduk(9);
	double sv = 0.001;
	
	double eu = get_loss(s0, sq, uk);
	cout << eu << endl;
	cout << "" << endl;
	for (int i = 0; i < 9; i++) {
		vec uk2 = uk;
		uk2(i) += sv;
		double eu2 = get_loss(s0, sq, uk2);
		deuduk(i) = (eu2 - eu) / sv;
	}
	return deuduk;
}

Mat<double> des_du(vec s0, vec uk) {
	Mat<double> desdu(5, 9);
	double sv = 1.0;

	vec es(5);
	es.zeros();
	for (int i = 0; i < 9; i++) {
		vec u = uk;
		u(i) += sv;
		vec sq = FK_3b(s0, u);
		vec uq = IK_3b(sq, uk, s0);
		es = FK_3b(s0, uq) - sq;
		desdu.col(i) = es;
	}
	return desdu;
}

double get_loss(vec s0, vec sq, vec uk) {
	vec uq = IK_3b(sq, uk, s0);
	//vec e = sq - FK_3b(s0, uq);
	Unicycle uc;
	uc.set_state(s0);
	uc.update_3b(uq);
	Mat<double> J = uc.jacobian(uq(0), uq(1), uq(2),
								uq(3), uq(4), uq(5), 
								uq(6), uq(7), uq(8));
	vec e = sq - J*(uq - uk) - FK_3b(s0, uk);
	vec uqr = uq - pinv(J) * e;
	vec eu = uqr - uk;

	double loss = dot(eu, eu.t());
	return loss;
}

/*
Mat<double> uniform_sample(vec s0, double da, double db, double dt) {
	int num = pow((1 + (A_MAX - A_MIN) / da), 3) * pow((1 + (B_MAX - B_MIN) / db), 3) * pow(1 + T_MAX / dt, 3);
}
*/


void kernel_state_test(vec s0, vec uk, bool thrb) {
	vec sk(5);
	if (!thrb) {
		sk = FK(s0, uk);
	}
	else {
		sk = FK_3b(s0, uk);
	}
	
	double itv = 0.01;
	double bound = 0.5;
	int num = pow(2.0 * bound / itv, 3);
	Mat<double> sq(5, num);
	Mat<double> uqs;
	if (thrb) {
		uqs = Mat<double>(9, num);
	}
	else {
		uqs = Mat<double>(3, num);
	}
	vec errors(num);
	errors.zeros();

	int idx = 0;
	for (int i = 0; i < 2.0 * bound / itv; i++) {
		for (int j = 0; j < 2.0 * bound / itv; j++) {
			for (int k = 0; k < 2.0 * bound / itv; k++) {
				// velocity and angular velocity
				sq(3, idx) = sk(3);
				sq(4, idx) = sk(4);
				// x, y, theta
				sq(0, idx) = sk(0) - bound + double(i) * itv;
				sq(1, idx) = sk(1) - bound + double(j) * itv;
				sq(2, idx) = sk(2) - bound + double(k) * itv;

				vec s(5);
				if (!thrb) {
					vec uq = IK(sq.col(idx), uk, s0);
					uqs.col(idx) = uq;
					s = FK(s0, uq);
				}
				else {
					vec uq = IK_3b(sq.col(idx), uk, s0);
					uqs.col(idx) = uq;
					s = FK_3b(s0, uq);
				}
				
				errors(idx) = euclidean_distance(s, sq.col(idx));
				idx++;
			}
		}
	}

	// Write Data
	ofstream myfile1, myfile2;
	myfile1.open("kernel_state.data");
	myfile2.open("kernel_u.data");
	myfile2 << "#uk is ";
	for (int i = 0; i < 9; i++) {
		myfile2 << uk(i) << " ";
	}
	myfile2 << "\n";

	for (int i = 0; i < num; i++) {
		for (int j = 0; j < 5; j++) {
			myfile1 << sq(j, i) << " ";
		}
		myfile1 << errors(i) << "\n";
		for (int j = 0; j < uqs.n_rows; j++) {
			myfile2 << uqs(j, i) << " ";
		}
		myfile2 << "\n";
	}
	myfile1.close();
	myfile2.close();
}

