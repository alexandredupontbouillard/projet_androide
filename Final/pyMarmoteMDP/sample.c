#include <stdio.h>
#include <math.h>
#include "law.h"

#define MAX_RANDOM_VALUE          0x7FFFFFFF

/* #ifndef ALPHA */
/* extern long random(void); */
/* #endif */

/****************************************************************/
/* Drawing of random numbers according to several classical	*/
/* distributions.						*/
/****************************************************************/
double Erlang(double N, double Mean)
/* N is kept double to avoid conversions */
{
  register int		i;
  register double	U;

  U = U_0_1();
	
  for (i=1; i<N; i++)	/* even if N not integer, this should	*/
			/* be the correct number of iterations	*/
    U *= U_0_1();

  return -Mean*log(U)/N;
}

double HyperExpo(double *Parameters)
{
  register double	U;
  register int    	i;

  U = U_0_1();

  for (i=0; (i<Parameters[0]) && (U>Parameters[i+1]); i+=2 )
    U -=Parameters[i+1];

  return Exponential( Parameters[i+2] );
}


/* Note: the generation of Weibull sambles would be faster if log(a)	*/
/* and 1/nu were stored instead of a and nu.				*/
double	Weibull(double a, double nu)
{
  register double	U;

  U = U_0_1();
  
  return pow(-log(U)/log(a),1.0/nu);
}

/* Generation of an exponential using the Polar Method (Knuth v2 p. 117) */
double Normal(double Mean, double StdDev )
{
  register double	v1;
  register double	v2;
  register double	S;

  /* Generate uniform on the unit circle */
  do { 
    v1 = 2 * U_0_1() - 1;
    v2 = 2 * U_0_1() - 1;
   } while ( 1 < ( S = v1*v1+v2*v2 ) );

  return Mean + StdDev * v1 * sqrt( -2.0 * log(S)/S );
}

/* Lognormal = exp(Normal). Adjusting the mean and var requires 	*/
/* computing parameters of the Normal. Would be faster if pre-computed.	*/
double LogNormal(double Mean, double StdDev )
{
  register double	m;
  register double	s2;

  s2 = log( 1.0 + (StdDev / Mean)*(StdDev / Mean) );
  m = log( Mean ) - s2 / 2;

  return exp( Normal( m, sqrt(s2) ) );
}

double Pareto( double a, double Alpha )
{
  return a * ( pow( U_0_1(), -1.0/Alpha ) - 1.0 );
}

double Sample(Law_Desc	The_Law)
{
  double	The_Val;
  static boolean	Sample_Error = FALSE;

  switch( The_Law->Name ) {
  case 'G':
    The_Val = 1.0 + Geometric( The_Law->Parameters[0] );
    break;
  case 'H':
    The_Val = HyperExpo( The_Law->Parameters );
    break;
  case 'I':
    The_Val = INFINITE_DURATION;
    break;
  case 'L':
    The_Val = LogNormal( The_Law->Parameters[0],
			 The_Law->Parameters[1] );
    break;
  case 'l':
    The_Val = LogNormal( The_Law->Parameters[0],
			 The_Law->Parameters[1] )
      + LogNormal( The_Law->Parameters[0],
		   The_Law->Parameters[1] );
    break;
  case 'N':
    The_Val = Normal( The_Law->Parameters[0],
		      The_Law->Parameters[1] );
    break;
  case 'p':
    The_Val = Pareto( The_Law->Parameters[0],
		      The_Law->Parameters[1] );
    break;
  case 'R':
    The_Val = Erlang( The_Law->Parameters[0],
		      The_Law->Parameters[1] );
    break;
  case 'W':
    The_Val = Weibull( The_Law->Parameters[0],
		       The_Law->Parameters[1] );
    break;
  default:
    if ( !Sample_Error ) {
      Sample_Error = TRUE;
      fprintf( stderr,
	       "Error: don't know hot to draw sample from distribution `%c`.",
	       The_Law->Name );
      fprintf( stderr,
	       " 0.0 assumed.\n");
    }
    The_Val = 0.0;
  }

  return The_Val;
}

void draw_iid_table(double	*The_Val,
		    int		The_Size,
		    Law_Desc	The_Law )
{	
  int	i;

  switch( The_Law->Name ) {
  case 'E':
    for (i=0; i<The_Size; i++)
      The_Val[i] = Exponential( The_Law->Parameters[0] );
    break;
  case 'g':
    for (i=0; i<The_Size; i++)
      The_Val[i] = Geometric( The_Law->Parameters[0] );
    break;
  case 'G':
    for (i=0; i<The_Size; i++)
      The_Val[i] = 1.0 + Geometric( The_Law->Parameters[0] );
    break;
  case 'H':
    for (i=0; i<The_Size; i++)
      The_Val[i] = HyperExpo( The_Law->Parameters );
    break;
  case 'I':
    for (i=0; i<The_Size; i++)
      The_Val[i] = INFINITE_DURATION;
    break;
  case 'L':
    for (i=0; i<The_Size; i++)
      The_Val[i] = LogNormal( The_Law->Parameters[0],
			      The_Law->Parameters[1] );
    break;
  case 'N':
    for (i=0; i<The_Size; i++)
      The_Val[i] = Normal( The_Law->Parameters[0],
			   The_Law->Parameters[1] );
    break;
  case 'p':
    for (i=0; i<The_Size; i++)
      The_Val[i] = Pareto( The_Law->Parameters[0],
			   The_Law->Parameters[1] );
    break;
  case 'R':
    for (i=0; i<The_Size; i++)
      The_Val[i] = Erlang( The_Law->Parameters[0],
			   The_Law->Parameters[1] );
    break;
  case 'W':
    for (i=0; i<The_Size; i++)
      The_Val[i] = Weibull( The_Law->Parameters[0],
			    The_Law->Parameters[1] );
  default:
    printf("Warning: unsupported law: %c. 0 values assumed.\n",
	   The_Law->Name );
    for (i=0; i<The_Size; i++)
      The_Val[i] = 0.0;		
  }
}
