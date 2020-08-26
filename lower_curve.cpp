/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#include "Curvebase.h"
#include "lower_curve.h"
#include <cmath>
lower_curve::lower_curve() {};
lower_curve::lower_curve(double p_min_in, double p_max_in, double d_in)
{
	pmin = p_min_in;
	pmax = p_max_in;
	length = integrate(pmin, pmax, tol);
	d = d_in;
};
double lower_curve::xp(double p) { return p; };
double lower_curve::dxp(double p) { return 1; };

double lower_curve::yp(double p)
{
	if (p < -3) {
		return 0.5 * 1 / (1 + exp(-3 * (p + 6)));
	}
	else {
		return 0.5 * 1 / (1 + exp(3 * (p)));
	}
};
double lower_curve::dyp(double p)
{
	if (p < -3) {
		return 3.0 / 2 * exp(-3 * (p)-18.0) / ((1 + exp(-3 * (p)-18)) * (1 + exp(-3 * (p)-18)));
	}
	else {
		return -3.0 / 2 * exp(3 * (p)) / ((1 + exp(3 * (p))) * (1 + exp(3 * (p))));
	}
};
double lower_curve::df(double p)
{
	return sqrt(dyp(p) * dyp(p) + dxp(p) * dxp(p));
};