#include "SPoint.h"
SPoint::SPoint()
{
	x = 0;
	y = 0;
	z = 0;
}
SPoint::~SPoint()
{
}
SPoint::SPoint(double x_, double y_, double z_)
{
	x = x_;
	y = y_;
	z = z_;
}