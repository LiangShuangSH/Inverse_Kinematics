#include "Sampler.h"

double acc_angle(vec s0, vec uk);
bool detect_singularity(vec uk);

void sampling(vec s0, double itv) {
	// Open file
	ofstream myfile_uk;
	ofstream myfile_sk;
	myfile_uk.open("sample_uk.data");
	myfile_sk.open("sample_sk.data");

	// Num of data
	int a_num = (A_MAX - A_MIN) / itv + 1;
	int b_num = (B_MAX - B_MIN) / itv + 1;
	int t_num = (T_MAX - T_MIN) / itv + 1;
	unsigned int data_num = pow(a_num * b_num * t_num, 3);

	// Uk initialization
	vec uk(9);

	// Sk initialization
	vec sk(5);

	for (unsigned int i = 0; i < data_num; i++) {
		int sample_idx[9] = {0};
		int num = i;
		int location = 0;

		while (num != 0) {
			int remainder;
			if (location % 3 == 0) {
				remainder = num % a_num;
				num = num / a_num;
			}
			else if (location % 3 == 1) {
				remainder = num % b_num;
				num = num / b_num;
			}
			else {
				remainder = num % t_num;
				num = num / t_num;
			}
			sample_idx[location] = remainder;
			location++;
		}

		// Translate to uk
		uk[0] = A_MIN + itv * sample_idx[0];
		uk[1] = B_MIN + itv * sample_idx[1];
		uk[2] = T_MIN + itv * sample_idx[2];
		uk[3] = A_MIN + itv * sample_idx[3];
		uk[4] = B_MIN + itv * sample_idx[4];
		uk[5] = T_MIN + itv * sample_idx[5];
		uk[6] = A_MIN + itv * sample_idx[6];
		uk[7] = B_MIN + itv * sample_idx[7];
		uk[8] = T_MIN + itv * sample_idx[8];


		// Write uk to file
		for (int j = 0; j < 9; j++) {
			myfile_uk << uk[j] << " ";
		}
		myfile_uk << "\n";

		// Write sk o file
		for (int j = 0; j < 5; j++) {
			myfile_sk << sk[j] << " ";
		}
		myfile_sk << "\n";

		// Show progress
		cout << int(float(i) / float(data_num) * 100.0) << "%\r";
		cout.flush();
	}

	myfile_uk.close();
	myfile_sk.close();
}

void uniform_sample_uspace(vec s0, Mat<double> &uks, Mat<double> &sks, double itv) {
	int a_num = (A_MAX - A_MIN) / itv + 1;
	int b_num = (B_MAX - B_MIN) / itv + 1;
	int t_num = (T_MAX - T_MIN) / itv + 1;
	unsigned int num = pow(a_num * b_num * t_num, 3);

	// Initializing
	uks = Mat<double>(9, num);
	sks = Mat<double>(5, num);

	vec uk(9);
	uk(0) = A_MIN; uk(1) = B_MIN; uk(2) = T_MIN;
	uk(3) = A_MIN; uk(4) = B_MIN; uk(5) = T_MIN;
	uk(6) = A_MIN; uk(7) = B_MIN; uk(8) = T_MIN;
	uks.col(0) = uk;
	vec sk(5);
	sk = FK_3b(s0, uk);
	sks.col(0) = sk;

	vec u_prev(9);
	u_prev = uk;

	for (unsigned int i = 1; i < num; i++) {
		uk = u_prev;
		uk(0) += itv;
		// Check carry
		uk(1) = (uk(0) > B_MAX) ? (uk(1) + itv) : uk(1);
		uk(2) = (uk(1) > T_MAX) ? (uk(2) + itv) : uk(2);
		uk(3) = (uk(2) > A_MAX) ? (uk(3) + itv) : uk(3);
		uk(4) = (uk(3) > B_MAX) ? (uk(4) + itv) : uk(4);
		uk(5) = (uk(4) > T_MAX) ? (uk(5) + itv) : uk(5);
		uk(6) = (uk(5) > A_MAX) ? (uk(6) + itv) : uk(6);
		uk(7) = (uk(6) > B_MAX) ? (uk(7) + itv) : uk(7);
		uk(8) = (uk(7) > T_MAX) ? (uk(8) + itv) : uk(8);
		for (int j = 0; j < 9; j++) {
			if (j % 3 == 0 && uk(j) > A_MAX) {
				uk(j) = A_MIN;
			}
			else if (j % 3 == 1 && uk(j) > B_MAX) {
				uk(j) = B_MIN;
			}
			else if (j % 3 == 2 && uk(j) > T_MAX) {
				uk(j) = T_MIN;
			}
		}

		u_prev = uk;

		uks.col(i) = uk;
		sk = FK_3b(s0, uk);
		sks.col(i) = sk;
	}
}

void write_data(Mat<double> uks, Mat<double> sks, bool overwrite) {
	ofstream myfile_uk;
	ofstream myfile_sk;
	if (overwrite) {
		myfile_uk.open("uniform_uk.data");
		myfile_sk.open("uniform_sk.data");
	}
	else {
		myfile_uk.open("uniform_uk.data", fstream::app);
		myfile_sk.open("uniform_sk.data", fstream::app);
	}

	int num_kernels = uks.n_cols;
	int u_dim = 9;
	int s_dim = 5;

	vec uk(u_dim);
	vec sk(s_dim);

	for (int i = 0; i < num_kernels; i++) {
		uk = uks.col(i);
		sk = sks.col(i);
		for (int j = 0; j < u_dim; j++) {
			myfile_uk << uk[j] << " ";
		}
		myfile_uk << "\n";
		for (int j = 0; j < s_dim; j++) {
			myfile_sk << sk[j] << " ";
		}
		myfile_sk << "\n";
	}
	myfile_uk.close();
	myfile_sk.close();
}

double acc_angle(vec s0, vec uk) {
	double theta = 0.0;
	int bands = uk.n_elem / 3;
	double* w = new double[bands+1];
	w[0] = s0[4];
	for (int i = 0; i < bands; i++) {
		double b = uk[3 * i + 1];
		double t = uk[3 * i + 2];
		w[i + 1] = w[i] + b*t;
		theta += w[i] * t + 0.5 * b * pow(t, 2);
	}
	
	w = NULL;
	delete w;

	return theta;
}

bool detect_singularity(vec uk) {
	bool sg = false;
	if (uk(2) + uk(5) + uk(8) == 0.0) {
		sg = true;
	}
	return sg;
}