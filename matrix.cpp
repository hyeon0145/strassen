#include "matrix.h"

int threshold = 64;

Matrix::Matrix(int number_of_rows, int number_of_columns) : elements(number_of_rows, Row(number_of_columns, 0)), number_of_rows(number_of_rows), number_of_columns(number_of_columns) {

}

Matrix::~Matrix() {

}

int Matrix::get(int row, int column) {
	assert(this->position_is_in_range(row, column));

	return this->elements[row][column];
}

void Matrix::set(int row, int column, int value) {
	assert(this->position_is_in_range(row, column));

	this->elements[row][column] = value;
}

Matrix Matrix::expand_row(int count) {
	Matrix result = Matrix(this->number_of_rows + count, this->number_of_columns);
	std::copy(this->elements.begin(), this->elements.end(), result.elements.begin());

	return result;
}

Matrix Matrix::expand_column(int count) {
	Matrix result = Matrix(this->number_of_rows, this->number_of_columns + count);
	for (int row = 0; row < this->number_of_rows; row += 1) {
		std::copy(this->elements[row].begin(), this->elements[row].end(), result.elements[row].begin());
	}

	return result;
}

Matrix Matrix::contract_row(int count) {
	Matrix result = Matrix(this->number_of_rows - count, this->number_of_columns);
	std::copy(this->elements.begin(), this->elements.end() - count, result.elements.begin());

	return result;
}

Matrix Matrix::contract_column(int count) {
	Matrix result = Matrix(this->number_of_rows, this->number_of_columns - count);
	for (int row = 0; row < this->number_of_rows; row += 1) {
		std::copy(this->elements[row].begin(), this->elements[row].end() - count, result.elements[row].begin());
	}

	return result;
}

Matrix Matrix::add(const Matrix &operand) {
	assert(this->has_same_size_with(operand));

	Matrix result(this->number_of_rows, this->number_of_rows);
	for (int row = 0; row < this->number_of_rows; row += 1) {
		for (int column = 0; column < this->number_of_columns; column += 1) {
			result.elements[row][column] = this->elements[row][column] + operand.elements[row][column];
		}
	}

	return result;
}

Matrix Matrix::subtract(const Matrix &operand) {
	assert(this->has_same_size_with(operand));

	Matrix result(this->number_of_rows, this->number_of_rows);
	for (int row = 0; row < this->number_of_rows; row += 1) {
		for (int column = 0; column < this->number_of_columns; column += 1) {
			result.elements[row][column] = this->elements[row][column] - operand.elements[row][column];
		}
	}

	return result;
}

Matrix Matrix::multiply(const Matrix &operand) {
	assert(this->is_able_to_multiply_with(operand));

	Matrix result(this->number_of_rows, operand.number_of_columns);
	for (int row = 0; row < this->number_of_rows; row += 1) {
		for (int column = 0; column < operand.number_of_columns; column += 1) {
			for (int n = 0; n < this->number_of_columns; n += 1) {
				result.elements[row][column] += this->elements[row][n] * operand.elements[n][column];
			}
		}
	}

	return result;
}

std::vector<Matrix> Matrix::partition() const {
	int partitioned_number_of_rows = this->number_of_rows / 2;
	int partitioned_number_of_columns = this->number_of_columns / 2;

	std::vector<Matrix> result;
	for (int i = 0; i < 4; i += 1) result.push_back(Matrix(partitioned_number_of_rows, partitioned_number_of_columns));

	for (int row = 0; row < partitioned_number_of_rows; row += 1) {
		std::copy(this->elements[row].begin(), this->elements[row].begin() + partitioned_number_of_columns, result[0].elements[row].begin());
		std::copy(this->elements[row].begin() + partitioned_number_of_columns, this->elements[row].begin() + (2 * partitioned_number_of_columns), result[1].elements[row].begin());
		std::copy(this->elements[partitioned_number_of_rows + row].begin(), this->elements[partitioned_number_of_rows + row].begin() + partitioned_number_of_columns, result[2].elements[row].begin());
		std::copy(this->elements[partitioned_number_of_rows + row].begin() + partitioned_number_of_columns, this->elements[partitioned_number_of_rows + row].begin() + (2 * partitioned_number_of_columns), result[3].elements[row].begin());
	}

	return result;
}

Matrix Matrix::multiply_by_strassen(const Matrix &operand) {
	assert(this->has_same_size_with(operand) && this->is_square_matrix());

	int size = this->number_of_rows; // size is equal to number_of_rows and number_of_columns of both lhs and rhs
	if (size <= threshold) {
		return this->multiply(operand);
	}

	Matrix lhs = *this, rhs = operand;

	bool is_expanded = false;
	if ((size % 2) == 1) {
		lhs = lhs.expand_row().expand_column();
		rhs = rhs.expand_row().expand_column();
		size += 1;

		is_expanded = true;
	}

	std::vector<Matrix> partitioned_lhs = lhs.partition();
	std::vector<Matrix> partitioned_rhs = rhs.partition();

	std::vector<Matrix> factors = std::vector<Matrix>(7);
	factors[0] = (partitioned_lhs[0] + partitioned_lhs[3]).multiply_by_strassen(partitioned_rhs[0] + partitioned_rhs[3]);
	factors[1] = (partitioned_lhs[2] + partitioned_lhs[3]).multiply_by_strassen(partitioned_rhs[0]);
	factors[2] = partitioned_lhs[0].multiply_by_strassen(partitioned_rhs[1] - partitioned_rhs[3]);
	factors[3] = partitioned_lhs[3].multiply_by_strassen(partitioned_rhs[2] - partitioned_rhs[0]);
	factors[4] = (partitioned_lhs[0] + partitioned_lhs[1]).multiply_by_strassen(partitioned_rhs[3]);
	factors[5] = (partitioned_lhs[2] - partitioned_lhs[0]).multiply_by_strassen(partitioned_rhs[0] + partitioned_rhs[1]);
	factors[6] = (partitioned_lhs[1] - partitioned_lhs[3]).multiply_by_strassen(partitioned_rhs[2] + partitioned_rhs[3]);

	std::vector<Matrix> elements = std::vector<Matrix>(4);
	elements[0] = factors[0] + factors[3] - factors[4] + factors[6];
	elements[1] = factors[2] + factors[4];
	elements[2] = factors[1] + factors[3];
	elements[3] = factors[0] + factors[2] - factors[1] + factors[5];
	assert((elements[0].number_of_rows == (size / 2)) && (elements[0].number_of_columns == (size / 2)));

	Matrix result(size, size);
	int size_of_element = (size / 2);
	for (int row = 0; row < size_of_element; row += 1) {
		std::copy(elements[0].elements[row].begin(), elements[0].elements[row].end(), result.elements[row].begin());
		std::copy(elements[1].elements[row].begin(), elements[1].elements[row].end(), result.elements[row].begin() + size_of_element);
		std::copy(elements[2].elements[row].begin(), elements[2].elements[row].end(), result.elements[size_of_element + row].begin());
		std::copy(elements[3].elements[row].begin(), elements[3].elements[row].end(), result.elements[size_of_element + row].begin() + size_of_element);
	}


	if (is_expanded) {
		result = result.contract_row().contract_column();
	}

	return result;
}

bool Matrix::is_equal_to(const Matrix &operand) {
	if ((this->number_of_rows != operand.number_of_rows)
		|| (this->number_of_columns != operand.number_of_columns)) return false;

	for (int row = 0; row < this->number_of_rows; row += 1) {
		for (int column = 0; column < this->number_of_columns; column += 1) {
			if (this->elements[row][column] != operand.elements[row][column]) return false;
		}
	}

	return true;
}

bool Matrix::row_is_in_range(int row) {
	return (0 <= row) && (row < this->number_of_rows);
}

bool Matrix::column_is_in_range(int column) {
	return (0 <= column) && (column < this->number_of_rows);
}

bool Matrix::position_is_in_range(int row, int column) {
	return this->row_is_in_range(row) && this->column_is_in_range(column);
}

bool Matrix::has_same_size_with(const Matrix &operand) {
	return (this->number_of_rows == operand.number_of_rows)
			&& (this->number_of_columns == operand.number_of_columns);
}

bool Matrix::is_square_matrix() {
	return (this->number_of_rows == this->number_of_columns);
}

bool Matrix::is_able_to_multiply_with(const Matrix &operand) {
	return (this->number_of_columns == operand.number_of_rows);
}

Matrix Matrix::operator +(const Matrix &operand) {
	return this->add(operand);
}

Matrix Matrix::operator -(const Matrix &operand) {
	return this->subtract(operand);
}

bool Matrix::operator ==(const Matrix &operand) {
	return this->is_equal_to(operand);
}

std::string Matrix::toString() {
	std::stringstream stream;

	for (int row = 0; row < this->number_of_rows; row += 1) {
		if (row != 0) stream << '\n';

		stream << "| ";
		for (int column = 0; column < this->number_of_columns; column += 1) {
			if (column != 0) stream << '\t';
			stream << this->elements[row][column];
		}
		stream << " |";
	}

	return stream.str();
}
