/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#include <cmath>
#include <list>
#include <iostream>
#include "Curvebase.h"
using namespace std;
Curvebase::Curvebase() {};
Curvebase::Curvebase(const Curvebase& object)
{
	pmin = object.pmin, pmax = object.pmax;
	length = object.length, d = object.d;
};
Curvebase::~Curvebase() {};
double Curvebase::x(double s)
{
	double p0 = 0, error = 1, toln = 1e-4, p, f;
	int iter = 0, iter_max = 500;
	while (error > toln&& iter < iter_max)
	{
		f = integrate(pmin, p0, tol) - s * length;
		p = p0 - f / df(p0);
		error = fabs(p - p0);
		p0 = p;
		iter += 1;
	}
	if (iter == iter_max)
	{
		cout << "Solution not yet converged to desired"
			" accuracy within 500 iterations." << endl;
	}
	return xp(p);
};
double Curvebase::y(double s)
{
	double p0 = 0, error = 1, toln = 1e-4, p, f;
	int iter = 0, iter_max = 500;
	while (error > toln&& iter < iter_max)
	{
		f = integrate(pmin, p0, tol) - s * length;
		p = p0 - f / df(p0);
		error = fabs(p - p0);
		p0 = p;
		iter += 1;
	}
	if (iter == iter_max)
	{
		cout << "Solution not yet converged to desired"
			" accuracy within 500 iterations." << endl;
	}
	return yp(p);
};
//declaration of function
double Curvebase::I_1(double pmin, double pmax)
{
	//definition of function, integral using Simpson's rule
	double I1 = ((pmax - pmin) / 6) * (df(pmin) + 4 * df((pmin + pmax) / 2) + df(pmax));
	return I1;
}
//declaration of function
double Curvebase::I_2(double pmin, double pmax)
{
	//definition of variable, midpoint
	double gamma = (pmin + pmax) * 0.5;
	//definition of function, helping calculate error
	double I2 = I_1(pmin, gamma) + I_1(gamma, pmax);
	return I2;
}
//declaration of function, adaptive Simpson's Integration
double Curvebase::integrate(double pmin, double pmax, double tol)
{
	//compute I1 integral between values pmin and pmax
	double I1 = I_1(pmin, pmax);
	//compute I2 integral for estimating error
	double I2 = I_2(pmin, pmax);
	//calculate error estimation
	double errest = abs(I1 - I2);

	if (errest < 15 * tol) //decide if error is too big
	{
		//end loop if error estimate satisfies tolerance
		return I2;
	}
	else
	{
		//keep looping if error estimate does not satisfy
		//condition, algorithm will subdivide interval of
		//integration in two and apply ASI recursively
		return integrate(pmin, (pmin + pmax) / 2, tol / 2)
			+ integrate((pmin + pmax) / 2, pmax, tol / 2);
	}
}

