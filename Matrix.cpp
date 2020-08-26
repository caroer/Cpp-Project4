#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include "Matrix.h"
#include <iostream>
#include <fstream>
using namespace std;

// Reference:
// https://codereview.stackexchange.com/questions/149669/c-operator-overloading-for-matrix-operations-follow-up


Matrix::Matrix(int n_rows, int n_cols)	// Constructor of class Matrix.
{
	rows = n_rows;
	cols = n_cols;
	Mat = new double* [rows];

	for (int i = 0; i < rows; i++)
	{
		Mat[i] = new double[cols]();		// For each row i every column value shall be set to 0.
	}
}

Matrix::Matrix(const Matrix& q)			// Copy constructor. 
{
	rows = q.rows;
	cols = q.cols;
	Mat = new double* [q.rows];
	for (int i = 0; i < q.rows; i++)
	{
		Mat[i] = new double[q.cols];
	}
	for (int i = 0; i < q.rows; i++) {
		for (int j = 0; j < q.cols; j++) {
			Mat[i][j] = q.Mat[i][j];
		}
	}
}


Matrix::~Matrix()		// Destructor. 
{
	for (int i = 0; i < rows; i++) {
		delete[] Mat[i];
	}
	delete[] Mat;
}


Matrix& Matrix::operator=(const Matrix& q)		// Overload operator. 
{
	if (this == &q)
		return(*this);
	delete[] Mat;
	this->rows = q.rows;
	this->cols = q.cols;
	this->Mat = new double* [rows];
	int ol = 0;

	for (int i = 0; i < rows; i++)
	{
		this->Mat[i] = new double[cols];
		ol++;
		for (int j = 0; j < cols; j++)
		{
			this->Mat[i][j] = q.Mat[i][j];
		}
	}
	return(*this);
}

Matrix Matrix::operator+(const Matrix& Q) const {
	Matrix temp(rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			temp.Mat[i][j] = Mat[i][j] + Q.Mat[i][j];
		}
	}
	return temp;
}

Matrix& Matrix::operator+=(const Matrix& Q)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			this->Mat[i][j] += Q.Mat[i][j];

		}
	}
	return (*this);
}


Matrix& Matrix::operator*=(const double num)
{
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			this->Mat[i][j] *= num;
		}
	}
	return *this;
}


Matrix Matrix::operator*(const double& num)const {
	Matrix temp(rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			temp.Mat[i][j] = Mat[i][j] * num;
		}
	}
	return temp;
}

void Matrix::fillMatrix(const int n_rows, const int n_cols, const Matrix& Q) {
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			cout << "Enter matrix value for row number " << i << " and column number " << j << ": ";
			cin >> Q.Mat[i][j];
		}
	}
}