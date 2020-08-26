/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 4: Caroline Eriksson
 * Date: December 2019
 * ===============================================================
 */
#include <iostream>
#include <fstream>
#include "lower_curve.h"
#include "line_v.h"
#include "line_h.h"
#include "Domain.h"
using namespace std;

Domain::Domain(lower_curve& curve_0_in, line_v& line_v1_in,
	line_h& line_h2_in, line_v& line_v3_in)
{
	// Allocating memory for: 
	// Grid points for boundary curve 0
	//curve_0 = lower_curve(curve_0_in);
	curve_0 = lower_curve(curve_0_in);
	x_0 = shared_ptr<double[]>(new double[n + 1]);
	y_0 = shared_ptr<double[]>(new double[n + 1]);
	// Grid points for boundary curve 1
	line_v1 = line_v(line_v1_in);
	x_1 = shared_ptr<double[]>(new double[m + 1]);
	y_1 = shared_ptr<double[]>(new double[m + 1]);
	// Grid points for boundary curve 2
	line_h2 = line_h(line_h2_in);
	x_2 = shared_ptr<double[]>(new double[n + 1]);
	y_2 = shared_ptr<double[]>(new double[n + 1]);
	// Grid points for boundary curve 3
	line_v3 = line_v(line_v3_in);
	x_3 = shared_ptr<double[]>(new double[m + 1]);
	y_3 = shared_ptr<double[]>(new double[m + 1]);
};

void Domain::grid_generation()
{
	double h1 = 1.0 / n; //Stepsize x-direction
	double h2 = 1.0 / m; //Stepsize y-direction
	// Where boundary curves will have n+1 grid points
	for (int i = 0; i < n + 1; i++)
	{
		// Generating grid points:
		// Grid points for boundary curve 0
		x_0[i] = curve_0.x(i * h1);
		y_0[i] = curve_0.y(i * h1);
		// Grid points for boundary curve 2
		x_2[i] = line_h2.x(i * h1);
		y_2[i] = line_h2.y(i * h1);
	}
	// Where boundary curves will have m+1 grid points
	for (int i = 0; i < m + 1; i++)
	{
		// Generating grid points:
		// Grid points for boundary curve 1
		x_1[i] = line_v1.x(i * h2);
		y_1[i] = line_v1.y(i * h2);
		// Grid points for boundary curve 3
		x_3[i] = line_v3.x(i * h2);
		y_3[i] = line_v3.y(i * h2);
	}
	int gridpoints_tot = 0;
	// Interpolating interior grid points using
	// The algebraic grid generation formula
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			x[i][j] = (1 - i * 1.0 / n) * x_3[j] + i * 1.0 / n * x_1[j]
				+ (1 - j * 1.0 / m) * x_0[i] + j * 1.0 / m * x_2[i]
				- (1 - i * 1.0 / n) * (1 - j * 1.0 / m) * (-10)
				- i * 1.0 / n * (1 - j * 1.0 / m) * (5)
				- (1 - i * 1.0 / n) * j * 1.0 / m * (-10)
				- i * 1.0 / n * j * 1.0 / m * (5);
			y[i][j] = (1 - i * 1.0 / n) * y_3[j] + i * 1.0 / n * y_1[j]
				+ (1 - j * 1.0 / m) * y_0[i] + j * 1.0 / m * y_2[i]
				- (1 - i * 1.0 / n) * (1 - j * 1.0 / m) * (0)
				- i * 1.0 / n * (1 - j * 1.0 / m) * (0)
				- (1 - i * 1.0 / n) * j * 1.0 / m * (3)
				- i * 1.0 / n * j * 1.0 / m * (3);
			gridpoints_tot += 1;
		}
	}
	cout << "The total amount of gridpoints is: " << gridpoints_tot << endl;
};

void Domain::writetofile() {    
	ofstream boundary_h;
	boundary_h.open("boundary_h.txt");
	for (int i = 0; i < n + 1; i++)
	{
		// Store boundary grid points in file:
		// Grid points for boundary curve 2
		boundary_h << x_2[i];
		boundary_h << "\t";
		boundary_h << y_2[i];
		boundary_h << "\t";
		// Grid points for boundary curve 0
		boundary_h << x_0[i];
		boundary_h << "\t";
		boundary_h << y_0[i];
		boundary_h << "\n";
	}
	boundary_h.close();
	ofstream boundary_v;
	boundary_v.open("boundary_v.txt");
	for (int i = 0; i < m + 1; i++) {
		// Store boundary grid points in file:
		// Grid points for boundary curve 3
		boundary_v << x_3[i];
		boundary_v << "\t";
		boundary_v << y_3[i];
		boundary_v << "\t";
		// Grid points for boundary curve 1
		boundary_v << x_1[i];
		boundary_v << "\t";
		boundary_v << y_1[i];
		boundary_v << "\n";
	}
	boundary_v.close();
	ofstream interior_x;
	interior_x.open("interior_x.txt");
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			// Store x - coordinates for interpolated
			// interior grid points in file:
			interior_x << x[i][j];
			interior_x << "\t";
		}
		interior_x << "\n";
	}
	interior_x.close();
	ofstream interior_y;
	interior_y.open("interior_y.txt");
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			// Store y - coordinates for interpolated
			// interior grid points in file:
			interior_y << y[i][j];
			interior_y << "\t";
		}
		interior_y << "\n";
	}
	interior_y.close();
};
