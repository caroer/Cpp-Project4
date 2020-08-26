/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#ifndef lowercurve_h
#define lowercurve_h
#include "Curvebase.h"
class lower_curve : public Curvebase
{
	friend class Domain;
public:
	lower_curve(double p_min_in, double p_max_in, double d_in);
	lower_curve();
protected:
	double df(double p);
	double xp(double p);
	double yp(double p);
	double dxp(double p);
	double dyp(double p);
};
#endif
