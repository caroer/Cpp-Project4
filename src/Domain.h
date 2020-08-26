/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 4: Caroline Eriksson
 * Date: December 2019
 * ===============================================================
 */
#include "line_h.h"
#include "line_v.h"
#include "lower_curve.h"
#include <memory>
using namespace std;

#ifndef DOMAIN_N
#define DOMAIN_N
class Domain
{
public:
	// Grid intervals
	const static int n = 49, m = 19;
	// Interpolated grid points
	double x[n + 1][m + 1], y[n + 1][m + 1];
	Domain(lower_curve& curve_0, line_v& line_v1,line_h& line_h2, line_v& line_v3); 
	// Function generating grid
	void grid_generation();
	int xsize() { return n; }
	int ysize() { return m; }
	//Point operator()(int i, int j);
	bool grid_valid() {
		if (n != 0 && m != 0) { return true; }
	}
	// Function writing the grid to a file
	void writetofile();
private:
	// Grid points for boundary curve 0
	lower_curve curve_0;
	shared_ptr<double[]> x_0; 
	shared_ptr<double[]> y_0;
	// Grid points for boundary curve 1
	line_v line_v1;
	shared_ptr<double[]> x_1;
	shared_ptr<double[]> y_1;
	// Grid points for boundary curve 2
	line_h line_h2;
	shared_ptr<double[]> x_2;
	shared_ptr<double[]> y_2;
	// Grid points for boundary curve 3
	line_v line_v3;
	shared_ptr<double[]> x_3;
	shared_ptr<double[]> y_3;
};
#endif
