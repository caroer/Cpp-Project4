/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#ifndef linev
#define linev
#include "Curvebase.h"
class line_v : public Curvebase
{
	friend class Domain;
public:
	line_v(double pmin_in, double pmax_in, double d_in);
	line_v();
protected:
	double df(double p);
	double xp(double p);
	double yp(double p);
	double dxp(double p);
	double dyp(double p);
};
#endif