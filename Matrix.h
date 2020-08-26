#include <iostream>
using namespace std;

#ifndef MATRIX_H
#define MATRIX_H

class Matrix {

public:
	int rows, cols;
	double** Mat;
	Matrix(int rows, int cols);
	Matrix(const Matrix&);
	~Matrix();
	Matrix& operator=(const Matrix&);
	Matrix operator+(const Matrix&) const;
	Matrix& operator+=(const Matrix&);
	Matrix& operator*=(const double);
	Matrix operator*(const double&) const;
	static void printMatrix(const Matrix& Q) {
		cout << "\nThe entered matrix is: " << endl;
		for (int i = 0; i < Q.rows; i++) {
			for (int j = 0; j < Q.cols; j++) {
				cout << Q.Mat[i][j] << ' ';
			}
			cout << ' ' << endl;
		}
	};
	void fillMatrix(const int, const int, const Matrix& Q);

};

#endif

