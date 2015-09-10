Levmar 2.5
http://www.ics.forth.gr/~lourakis/levmar/


In levmar.h: uncomment #undef HAVE_LAPACK

This means that Lapack is not available and only unbounded problems can be solved. Furthermore LU will be used for matrix inversion when computing the covariance; this might be unstable at times.

Comment out #define LM_SNGL_PREC when only using double precision data

In compiler.h:
Change
//#define LM_FINITE isfinite // other than MSVC, ICC, GCC, let's hope this will work
to
#define LM_FINITE(x) finite_function(x)
and add 
int finite_function(double x);

Create compiler.c
#include "toolbox.h"
int finite_function(double x)
{
	if ( IsInfinity(x) || IsNotANumber(x) ) return 0;
	return 1;
}

Comment out a #warning in misc_core that the compiler complains about

Add suffix fit_ to files and replace references to filenames in .c and .h files by trying to compile