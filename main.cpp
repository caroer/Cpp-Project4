/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 4: Caroline Eriksson
 * Date: December 2019
 * ===============================================================
 */
#include "lower_curve.h"
#include "line_v.h"
#include "line_h.h"
#include "Domain.h"
#include "Matrix.h"
#include "GFkt.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h>
using namespace std;

int main()
{
	// Values for the boundary sides
	double pmin0 = -10.0, pmax0 = 5.0, d0 = 10.0;
	double pmin1 = 0.0, pmax1 = 3.0, d1 = 5.0;
	double pmin2 = -10.0, pmax2 = 5.0, d2 = 3.0;
	double pmin3 = 0.0, pmax3 = 3.0, d3 = -10.0;

	// Creation of the boundary sides
	lower_curve curve_0 = lower_curve(pmin0, pmax0, d0);
	line_v line_v1 = line_v(pmin1, pmax1, d1);
	line_h line_h2 = line_h(pmin2, pmax2, d2);
	line_v line_v3 = line_v(pmin3, pmax3, d3);

	// Creation of the domain and corresponding grid
	shared_ptr<Domain> domain = shared_ptr<Domain>(new Domain(curve_0, line_v1, line_h2, line_v3));
	domain->grid_generation();
	
	// Creation of the "grid function object"
	GFkt gf = GFkt(domain);

	// Discretizing the function u
	for (int i = 0; i < domain->n + 1; i++) {
		for (int j = 0; j < domain->m + 1; j++) {
			gf.u.Mat[i][j] = sin(pow(((domain->x[i][j]) / 10), 2)) * cos((domain->x[i][j]) / 10) + domain->y[i][j];
		}
	}
	
	cout << gf.u.Mat[14][0] << " " << gf.u.Mat[12][0] << endl;
	// Creation of the first order derivatives

	GFkt d_x = gf.D0x();
	GFkt d_y = gf.D0y();
	
	// Creation of the second order derivatives
	GFkt d2_x = gf.D2x();
	GFkt d2_y = gf.D2y();

	// Write the grid, the derivatives and the laplacian to files
	domain->writetofile();
	d_x.writetofile_dx();
	d_y.writetofile_dy();
	GFkt::writetofile_laplace(d2_x, d2_y, d2_x.grid);

	return 0;
}