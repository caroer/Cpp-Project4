/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#include "Curvebase.h"
#include "line_v.h"
line_v::line_v() {};
line_v::line_v(double p_min_in, double p_max_in, double d_in)
{
	pmin = p_min_in;
	pmax = p_max_in;
	length = p_max_in - p_min_in;
	d = d_in;
};
double line_v::xp(double p) { return d; };
double line_v::dxp(double p) { return 0.0; };
double line_v::yp(double p) { return p; };
double line_v::dyp(double p) { return 1.0; };
double line_v::df(double p) { return 1.0; };
