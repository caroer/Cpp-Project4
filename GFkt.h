/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 4: Caroline Eriksson
 * Date: December 2019
 * ===============================================================
 */
#include "Domain.h"
#include "Matrix.h"
#include <stdexcept>
#include <memory>

#ifndef GFKT_H
#define GFKT_H

class GFkt {
public:
	Matrix u;
	shared_ptr<Domain> grid;
	// Standard operations
	GFkt(shared_ptr<Domain> grid_) :
		u(grid_->xsize() + 1, grid_->ysize() + 1),
		grid(grid_) {}
	GFkt(const GFkt& U) : u(U.u), grid(U.grid) {}
	GFkt& operator = (const GFkt & U);
	GFkt operator+(const GFkt& U) const;
	GFkt operator*(const GFkt& U) const;
	// First and second order differentiations 
	GFkt D0x() const;
	GFkt D0y() const;
	GFkt D2x() const;
	GFkt D2y() const;
	// Writing to files
	void writetofile_dx();
	void writetofile_dy();
	static void writetofile_laplace(GFkt& d2_x, GFkt& d2_y, shared_ptr<Domain> grid);
};

#endif
