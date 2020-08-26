/* ===============================================================
 * SF2565 Program construction in C++ for Scientific Computing
 * Project 3: Sophie Malmliden & Caroline Eriksson
 * Date: November 2019
 * ===============================================================
 */
#ifndef lineh
#define lineh
#include "Curvebase.h"
class line_h : public Curvebase
{
	friend class Domain;
public:
	line_h(double p_min_in, double p_max_in, double d_in);
	line_h();
protected:
	double df(double p);
	double xp(double p);
	double yp(double p);
	double dxp(double p);
	double dyp(double p);
};
#endif

