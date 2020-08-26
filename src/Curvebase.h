/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#ifndef CURVEBASE
#define CURVEBASE
#include <cmath>
class Curvebase
{
protected:
	double pmin;               // Minimal value for p
	double pmax;               // Maximal value for p
	double length;
	double d;
	double tol = pow(10, -2);
	virtual double xp(double p) = 0;
	virtual double yp(double p) = 0;
	virtual double dxp(double p) = 0;
	virtual double dyp(double p) = 0;
	// Arc length integral
	double integrate(double pmin, double pmax, double tol);
	virtual double df(double) = 0;
public:
	Curvebase();                   // Constructor
	~Curvebase();                  // Destructor
	Curvebase(const Curvebase& object); // Copy constructor
	double x(double s);            // Arc length parametrization
	double y(double s);            // Arc length parametrization
	double I_1(double, double);    // Integral used in ASI
	double I_2(double, double);    // Integral used in ASI
};
#endif
