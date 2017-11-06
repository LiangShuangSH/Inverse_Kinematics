#include "Experiment.h"

void Experiment() {
	int num = 1000;
	int idx = 0;

	vec s0(9);
	s0.zeros();

	Mat<double> uk1_record(9, num);
	Mat<double> uk_true_record(9, num);
	Mat<double> sq_record(5, num);
	Mat<double> uk2_record(9, num);
	Mat<double> sd_record(5, num);
	
	vec error_record(num);

	// Do main stuff
	while (idx < num) {
		vec uk1(9);
		vec uk_true(9);
		vec uk2(9);
		vec sq(5);
		vec s2(5);
		double error;

		uk1 = random_uk(9);
		// If uk1 results whirling, then pass
		if (check_conditions(s0, uk1)) {
			continue;
		}
		vec sk = FK_3b(s0, uk1);

		uk_true = uk1 + normalize_vector(random_uk(9), 9, 0.1);
		sq = FK_3b(s0, uk_true);
		double distance = euclidean_distance(sq, sk);
		
		uk2 = IK_3b(sq, uk1, s0);

		// If t1, t2, t3, anyone of them <0.0, cut it out
		if (uk2(2) < 0.0 || uk2(5) < 0.0 || uk2(8) < 0.0) {
			continue;
		}

		s2 = FK_3b(s0, uk2);

		vec sd = sq - sk;

		error = euclidean_distance(s2, sq);

		// Recording
		uk1_record.col(idx) = uk1;
		uk_true_record.col(idx) = uk_true;
		sq_record.col(idx) = sq;
		uk2_record.col(idx) = uk2;
		sd_record.col(idx) = sd;
		error_record(idx) = error;

		idx += 1;
		std::cout << int(float(idx) / float(num) * 100.0) << "%\r";
	}

	uvec idx_sort = sort_index(error_record, "descend");

	// Write file
	ofstream error_file, uk1_file, uk_true_file, sq_file, uk2_file, sd_file;
	error_file.open("./Experiment/error.data");
	uk1_file.open("./Experiment/uk1.data");
	uk_true_file.open("./Experiment/uk_true.data");
	sq_file.open("./Experiment/sq.data");
	uk2_file.open("./Experiment/uk2.data");
	sd_file.open("./Experiment/sd.data");

	for (int i = 0; i < idx_sort.size(); i++) {
		int index = idx_sort(i);

		//Write error
		error_file << error_record(index) << "\n";

		// Write sq
		for (int j = 0; j < 5; j++) {
			sq_file << sq_record(j, index) << " ";
		}
		sq_file << "\n";

		// Write sd
		for (int j = 0; j < 5; j++) {
			sd_file << sd_record(j, index) << " ";
		}
		sd_file << "\n";

		// Write u
		for (int j = 0; j < 9; j++) {
			uk1_file << uk1_record(j, index) << " ";
			uk_true_file << uk_true_record(j, index) << " ";
			uk2_file << uk2_record(j, index) << " ";
		}
		uk1_file << "\n";
		uk_true_file << "\n";
		uk2_file << "\n";
	}

	error_file.close();
	uk1_file.close();
	uk_true_file.close();
	sq_file.close();
	uk2_file.close();
	sd_file.close();

	return;
}

void Experiment2() {
	int num = 100;
	int idx = 0;

	vec s0(9);
	s0.zeros();

	Mat<double> uk1_record(9, num);
	Mat<double> sq_record(5, num);
	Mat<double> uk2_record(9, num);
	Mat<double> sd_record(5, num);

	vec error_record(num);

	// Do main stuff
	while (idx < num) {
		vec uk1(9);
		vec uk2(9);
		vec sq(5);
		vec s2(5);
		double error;

		uk1 = random_uk(9);
		// If uk1 results whirling, then pass
		if (check_conditions(s0, uk1)) {
			continue;
		}
		vec sk = FK_3b(s0, uk1);

		sq = sk + normalize_vector(random_s(), 5, 0.1);
		double distance = euclidean_distance(sq, sk);

		uk2 = IK_3b(sq, uk1, s0);

		// If t1, t2, t3, anyone of them <0.0, cut it out
		if (uk2(2) < 0.0 || uk2(5) < 0.0 || uk2(8) < 0.0) {
			continue;
		}

		s2 = FK_3b(s0, uk2);

		vec sd = sq - sk;

		error = euclidean_distance(s2, sq);

		// Recording
		uk1_record.col(idx) = uk1;
		sq_record.col(idx) = sq;
		uk2_record.col(idx) = uk2;
		sd_record.col(idx) = sd;
		error_record(idx) = error;

		idx += 1;
		std::cout << int(float(idx) / float(num) * 100.0) << "%\r";
	}

	uvec idx_sort = sort_index(error_record, "descend");

	// Write file
	ofstream error_file, uk1_file, sq_file, uk2_file, sd_file;
	error_file.open("./Experiment/error_2.data");
	uk1_file.open("./Experiment/uk1_2.data");
	sq_file.open("./Experiment/sq_2.data");
	uk2_file.open("./Experiment/uk2_2.data");
	sd_file.open("./Experiment/sd_2.data");

	for (int i = 0; i < idx_sort.size(); i++) {
		int index = idx_sort(i);
		//Write error
		error_file << error_record(index) << "\n";

		// Write sq
		for (int j = 0; j < 5; j++) {
			sq_file << sq_record(j, index) << " ";
		}
		sq_file << "\n";

		// Write sd
		for (int j = 0; j < 5; j++) {
			sd_file << sd_record(j, index) << " ";
		}
		sd_file << "\n";

		// Write u
		for (int j = 0; j < 9; j++) {
			uk1_file << uk1_record(j, index) << " ";
			uk2_file << uk2_record(j, index) << " ";
		}
		uk1_file << "\n";
		uk2_file << "\n";
	}

	error_file.close();
	uk1_file.close();
	sq_file.close();
	uk2_file.close();
	sd_file.close();

	return;
}


void sample(int num) {
	int idx = 0;

	ofstream random_data, s_data;
	random_data.open("./Data/random_sample.data", fstream::app);
	s_data.open("./Data/s.data", fstream::app);

	while (idx < num) {
		// Random s0
		vec s0(9);
		s0.zeros();

		s0(3) = (-V_MAX) + ((double)rand() / RAND_MAX) * (V_MAX + V_MAX);
		s0(4) = (-W_MAX) + ((double)rand() / RAND_MAX) * (W_MAX + W_MAX);

		// Random u
		vec u = random_uk(9);

		// Check conditions
		if (check_conditions(s0, u)) {
			continue;
		}

		vec s1 = FK_3b(s0, u);

		// Write sample data
		random_data << s0(3) << " " << s0(4) << " ";

		for (int i = 0; i < 9; i++) {
			random_data << u(i) << " ";
		}

		random_data << "\n";

		// Write s data
		for (int i = 0; i < 5; i++) {
			s_data << s1(i) << " ";
		}
		s_data << "\n";

		idx += 1;
		std::cout << int(float(idx) / float(num) * 100.0) << "%\r";
	}

	random_data.close();
	s_data.close();

	return;
}

bool check_conditions(vec s0, vec uk) {
	Unicycle uc;
	uc.set_state(s0);

	double angle = 0.0;

	double a1 = uk(0);
	double a2 = uk(3);
	double a3 = uk(6);

	double b1 = uk(1);
	double b2 = uk(4);
	double b3 = uk(7);

	double t1 = uk(2);
	double t2 = uk(5);
	double t3 = uk(8);

	double delta_t = 0.01;

	double t = 0.0;
	while (t < t1) {
		uc.update(a1, b1, delta_t);
		vec s = uc.get_state_vec();
		double z = s(2);
		angle += z;
		// Check whiriling
		if (abs(angle) > PII) {
			return true;
		}
		// Check velocity
		double v = s(3);
		if (abs(v) > V_MAX) {
			return true;
		}
		// Check angular velocity
		double w = s(4);
		if (abs(w) > W_MAX) {
			return true;
		}
		t += delta_t;
	}

	while (t < t1 + t2) {
		uc.update(a2, b2, delta_t);
		vec s = uc.get_state_vec();
		double z = s(2);
		angle += z;
		// Check whiriling
		if (abs(angle) > PII) {
			return true;
		}
		// Check velocity
		double v = s(3);
		if (abs(v) > V_MAX) {
			return true;
		}
		// Check angular velocity
		double w = s(4);
		if (abs(w) > W_MAX) {
			return true;
		}
		t += delta_t;
	}

	while (t < t1 + t2 + t3) {
		uc.update(a3, b3, delta_t);
		vec s = uc.get_state_vec();
		double z = s(2);
		angle += z;
		// Check whiriling
		if (abs(angle) > PII) {
			return true;
		}
		// Check velocity
		double v = s(3);
		if (abs(v) > V_MAX) {
			return true;
		}
		// Check angular velocity
		double w = s(4);
		if (abs(w) > W_MAX) {
			return true;
		}
		t += delta_t;
	}

	return false;
}
