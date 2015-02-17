#pragma once

#include <iostream>

#include <vector>
#include <string>
#include <sstream>
#include <cassert>

typedef std::vector<int> Row;

extern int threshold;

class Matrix {
private:
	std::vector<Row> elements;

	int number_of_rows;
	int number_of_columns;

public:
	Matrix(int number_of_rows = 0, int number_of_columns = 0);
	virtual ~Matrix();

	int get(int row, int column);
	void set(int row, int column, int value);

	Matrix expand_row(int count = 1);
	Matrix expand_column(int count = 1);

	Matrix contract_row(int count = 1);
	Matrix contract_column(int count = 1);

	Matrix add(const Matrix &operand);
	Matrix subtract(const Matrix &operand);
	Matrix multiply(const Matrix &operand);
private:
	std::vector<Matrix> partition() const;
public:
	Matrix multiply_by_strassen(const Matrix &operand);

	bool is_equal_to(const Matrix &operand);

	inline bool row_is_in_range(int row);
	inline bool column_is_in_range(int column);
	inline bool position_is_in_range(int row, int column);

	inline bool has_same_size_with(const Matrix &operand);
	inline bool is_square_matrix();

	inline bool is_able_to_multiply_with(const Matrix &operand);

	inline Matrix operator +(const Matrix &operand);
	inline Matrix operator -(const Matrix &operand);
	bool operator ==(const Matrix &operand);

	std::string toString();
};

