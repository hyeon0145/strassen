#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cassert>
#include <limits>

#include "matrix.h"

Matrix create_matrix(int n) {
	Matrix result = Matrix(n, n);
	for (int row = 0; row < n; row += 1) {
		for (int column = 0; column < n; column += 1) {
			result.set(row, column, std::rand() % 1000);
		}
	}

	return result;
}

#define MATRIX_SIZE 1000

int main() {
	std::srand(std::time(0));

	Matrix lhs = create_matrix(MATRIX_SIZE);
	Matrix rhs = create_matrix(MATRIX_SIZE);

	int efficient_threshold = 1;
	double minimum_elapsed_time = std::numeric_limits<double>::max();
	for (int t = 1; t <= 512; t += 1) {
		threshold = t;
		std::cout << "threshold = " << threshold << std::endl;

		clock_t start_time = std::clock();
		Matrix result = lhs.multiply_by_strassen(rhs);
		clock_t end_time = std::clock();

		double elapsed_time = (end_time - start_time) / (double)CLOCKS_PER_SEC;
		std::cout << "elapsed time: " << elapsed_time << "sec" << std::endl;

		if (elapsed_time < minimum_elapsed_time) {
			efficient_threshold = threshold;
			minimum_elapsed_time = elapsed_time;
		}
	}

	std::cout << "efficient threshold is " << efficient_threshold << std::endl;

	return 0;
}
