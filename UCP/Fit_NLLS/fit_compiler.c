#include "toolbox.h"
int finite_function(double x)
{
	if ( IsInfinity(x) || IsNotANumber(x) ) return 0;
	return 1;
}