#include "Kinematics.h"

vec FK(vec s0, vec u) {
	Unicycle uc;
	uc.set_state(s0);
	uc.update(u);
	vec s = uc.get_state_vec();
	return s;
}

vec FK_3b(vec s0, vec u) {
	Unicycle uc;
	uc.set_state(s0);
	uc.update_3b(u);
	vec s = uc.get_state_vec();
	return s;
}

vec IK(vec sq, vec uk, vec s0) {
	double a = uk(0);
	double b = uk(1);
	double t = uk(2);

	Unicycle uc;
	uc.set_state(s0);
	uc.update(a, b, t);
	Mat<double> jcb = uc.jacobian(a, b, t);
	vec sk = uc.get_state_vec();

	Mat<double> jcb_inv = pinv(jcb);
	vec uq = jcb_inv * (sq - sk) + uk;
	return uq;
}

vec IK_3b(vec sq, vec uk, vec s0) {
	double a1 = uk(0);
	double b1 = uk(1);
	double t1 = uk(2);
	double a2 = uk(3);
	double b2 = uk(4);
	double t2 = uk(5);
	double a3 = uk(6);
	double b3 = uk(7);
	double t3 = uk(8);

	Unicycle uc;
	uc.set_state(s0);
	uc.update_3b(uk);

	Mat<double> jcb = uc.jacobian(a1, b1, t1, a2, b2, t2, a3, b3, t3);
	vec sk = uc.get_state_vec();

	Mat<double> jcb_inv = pinv(jcb);
	
	/*
	Mat<double> subJcb1 = jcb.submat(0, 0, 4, 2);
	Mat<double> subJcb2 = jcb.submat(0, 3, 4, 5);
	Mat<double> subJcb3 = jcb.submat(0, 6, 4, 8);

	Mat<double> invJcb1 = pinv(subJcb1);
	Mat<double> invJcb2 = pinv(subJcb2);
	Mat<double> invJcb3 = pinv(subJcb3);

	Mat<double> jcb_inv = join_cols(join_cols(invJcb1, invJcb2), invJcb3);*/

	vec uq = jcb_inv * (sq - sk) + uk;
	return uq;
}

vec random_uk(int dim) {
	vec uk(dim);

	for (int i = 0; i < dim; i++) {
		if (i%3 == 0) {
			uk(i) = A_MIN + ((double)rand() / RAND_MAX) * (A_MAX - A_MIN);
		}
		else if (i%3 == 1) {
			uk(i) = B_MIN + ((double)rand() / RAND_MAX) * (B_MAX - B_MIN);
		}
		else if (i%3 == 2) {
			uk(i) = 0.0 + ((double)rand() / RAND_MAX) * (T_MAX - 0.0);
		}
	}

	return uk;
}

vec random_s() {
	vec s(5);

	for (int i = 0; i < 5; i++) {
		s(i) = -1.0 + ((double)rand() / RAND_MAX) * (1.0 - (-1.0));
	}

	return s;
}