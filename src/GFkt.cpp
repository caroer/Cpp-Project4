/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 4: Caroline Eriksson
 * Date: December 2019
 * ===============================================================
 */

#include "GFkt.h"
#include "Domain.h"
#include "Matrix.h"
#include <iostream>
#include <fstream>

GFkt& GFkt::operator=(const GFkt& U) {
	if (this == &U) {
		return(*this);
	}
	this->u = U.u;
	this->grid = U.grid;
	return (*this);
}

GFkt GFkt::operator+(const GFkt& U) const {
	if (grid == U.grid) { // defined on the same grid?
		GFkt tmp(grid);
		tmp.u = u + U.u; // Matrix::operator+()
		return tmp;
	}
	else throw std::runtime_error("Not defined on the same grid.");
}
GFkt GFkt::operator*(const GFkt& U) const {
	if (grid == U.grid) { // defined on the same grid?
		GFkt tmp(grid);
		for (int j = 0; j <= grid->ysize(); j++)
			for (int i = 0; i <= grid->xsize(); i++)
				tmp.u.Mat[i][j] = u.Mat[i][j] * U.u.Mat[i][j];
				
		return tmp;
	}
	else throw std::runtime_error("Not defined on the same grid.");
}/* First order differentiation with respect to x */
GFkt GFkt::D0x() const {
	GFkt tmp(grid);
	double v1, v2;
	double n = grid->n;
	double dx = 15 / n;	
	if (grid->grid_valid()) {
		//Generating derivative in tmp
		for (int j = 0; j <= grid->ysize(); j++) {
			for (int i = 0; i <= grid->xsize(); i++) {
				if (i == 0) {		// Ghost points occuring
					v1 = sin(pow(((grid->x[i][j] - dx) / 10), 2)) * cos((grid->x[i][j] - dx) / 10) + grid->y[i][j];
					tmp.u.Mat[i][j] = (u.Mat[i + 1][j] - v1) / (2 * dx);
					continue;
				}
				if (i == grid->xsize()) {	// Ghost points occuring
					v2 = sin(pow(((grid->x[i][j] + dx) / 10), 2)) * cos((grid->x[i][j] + dx) / 10) + grid->y[i][j];
					tmp.u.Mat[i][j] = (v2 - u.Mat[i - 1][j]) / (2 * dx);
					continue;
				}
				// Inside the boundaries where no ghost points occurs
				//if (x[i][j])
				tmp.u.Mat[i][j] = (u.Mat[i + 1][j] - u.Mat[i - 1][j]) / (2 * dx);
			}
		}
	}
	return tmp;
}

/* First order differentiation with respect to y */GFkt GFkt::D0y() const {
	GFkt tmp(grid);
	double v1, v2;
	double m = grid->m;
	double dy = 3 / m;		
	if (grid->grid_valid()) {
		// Generating derivative in tmp
		for (int i = 0; i <= grid->xsize(); i++) {
			for (int j = 0; j <= grid->ysize(); j++) {
				if (j == 0) {		// Ghost points occuring
					v1 = sin(pow(((grid->x[i][j]) / 10), 2)) * cos((grid->x[i][j]) / 10) + grid->y[i][j] - dy;
					tmp.u.Mat[i][j] = (u.Mat[i][j + 1] - v1) / (2 * dy);
					continue;
				}
				if (j == grid->ysize()) {	// Ghost points occuring
					v2 = sin(pow(((grid->x[i][j]) / 10), 2)) * cos((grid->x[i][j]) / 10) + grid->y[i][j] + dy;
					tmp.u.Mat[i][j] = (v2 - u.Mat[i][j - 1]) / (2 * dy);
					continue;
				}
				// Inside the boundaries where no ghost points occurs
				tmp.u.Mat[i][j] = (u.Mat[i][j + 1] - u.Mat[i][j - 1]) / (2 * dy);
			}
		}
	}	return tmp;}

/* Second order differentiation with respect to x */
GFkt GFkt::D2x() const {
	GFkt tmp(grid);
	double v1, v2;
	double n = grid->n;
	double dx = 15 / n;	
	if (grid->grid_valid()) {
		// Generating derivative in tmp
		for (int j = 0; j <= grid->ysize(); j++) {
			for (int i = 0; i <= grid->xsize(); i++) {
				if (i == 0) {		// Ghost points occuring
					v1 = sin(pow(((grid->x[i][j] - dx) / 10), 2)) * cos((grid->x[i][j] - dx) / 10) + grid->y[i][j];
					tmp.u.Mat[i][j] = (u.Mat[i + 1][j] - 2*u.Mat[i][j] + v1) / (pow(dx,2));
					continue;
				}
				if (i == grid->xsize()) {	// Ghost points occuring
					v2 = sin(pow(((grid->x[i][j] + dx) / 10), 2)) * cos((grid->x[i][j] + dx) / 10) + grid->y[i][j];
					tmp.u.Mat[i][j] = (v2 - 2 * u.Mat[i][j] + u.Mat[i - 1][j]) / (pow(dx, 2));
					continue;
				}
				// Inside the boundaries where no ghost points occurs
				tmp.u.Mat[i][j] = (u.Mat[i + 1][j] - 2 * u.Mat[i][j] + u.Mat[i - 1][j]) / (pow(dx, 2));
			}
		}
	}
	return tmp;
}

/* Second order differentiation with respect to y */
GFkt GFkt::D2y() const {
	GFkt tmp(grid);
	double v1, v2;
	double m = grid->m;
	double dy = 3 / m;
	if (grid->grid_valid()) {
		// Generating derivative in tmp
		for (int i = 0; i <= grid->xsize(); i++) {
			for (int j = 0; j <= grid->ysize(); j++) {
				if (j == 0) {		// Ghost points occuring
					v1 = sin(pow(((grid->x[i][j]) / 10), 2)) * cos((grid->x[i][j]) / 10) + grid->y[i][j] - dy;
					tmp.u.Mat[i][j] = (u.Mat[i][j + 1] - 2 * u.Mat[i][j] + v1) / (pow(dy, 2));
					continue;
				}
				if (j == grid->ysize()) {	// Ghost points occuring
					v2 = sin(pow(((grid->x[i][j]) / 10), 2)) * cos((grid->x[i][j]) / 10) + grid->y[i][j] + dy;
					tmp.u.Mat[i][j] = (v2 - 2 * u.Mat[i][j] + u.Mat[i][j - 1]) / (pow(dy, 2));
					continue;
				}
				// Inside the boundaries where no ghost points occurs
				tmp.u.Mat[i][j] = (u.Mat[i][j + 1] - 2 * u.Mat[i][j] + u.Mat[i][j - 1]) / (pow(dy, 2));
			}
		}
	}	return tmp;}

/*Writing dx to a file */
void GFkt::writetofile_dx() {
	ofstream differential_x;
	differential_x.open("differential_x.txt");
	for (int i = 0; i < grid->xsize() + 1; i++)
	{
		for (int j = 0; j < grid->ysize() + 1; j++)
		{
			differential_x << u.Mat[i][j];
			differential_x << "\t";
		}
		differential_x << "\n";
	}
	differential_x.close();
};

/*Writing dy to a file */
void GFkt::writetofile_dy() {
	ofstream differential_y;
	differential_y.open("differential_y.txt");
	for (int i = 0; i < grid->xsize() + 1; i++)
	{
		for (int j = 0; j < grid->ysize() + 1; j++)
		{
			differential_y << u.Mat[i][j];
			differential_y << "\t";
		}
		differential_y << "\n";
	}
	differential_y.close();
};

/*Computing and writing the laplacian to a file */
void GFkt::writetofile_laplace(GFkt& d2_x, GFkt& d2_y, shared_ptr<Domain> grid) {
	GFkt tmp(grid);
	if (grid->grid_valid()) {
		// Computation of the laplacian
		for (int i = 0; i <= grid->xsize(); i++) {
			for (int j = 0; j <= grid->ysize(); j++) {
				tmp.u.Mat[i][j] = d2_x.u.Mat[i][j] + d2_y.u.Mat[i][j];
			}
		}
	}
	// Writing the laplacian to a file
	ofstream laplace;
	laplace.open("laplace.txt");
	for (int i = 0; i <= grid->xsize(); i++)
	{
		for (int j = 0; j <= grid->ysize(); j++)
		{
			laplace << tmp.u.Mat[i][j];
			laplace << "\t";
		}
		laplace << "\n";
	}
	laplace.close();
};

