/* Marmote is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Marmote is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Marmote. If not, see <http://www.gnu.org/licenses/>.

Copyright 2015 Alain Jean-Marie, Jean-Michel Fourneau, Jean-Marc Vincent, Issam Rabhi */

#include "geometricDistribution.h"
#include <stdlib.h>
#include <sstream>

/**
 * Constructor for a geometric distribution
 * The mean is calculated at creation.
 * The value p = 1 is admitted, in which case the law is a Dirac at infinity
 * and has no moments.
 * @author	Alain Jean-Marie
 * @param p the probability of being non 0
 * @return	an object of type geometricDistribution
 */
geometricDistribution::geometricDistribution( double p )
{
  _name = "geometricDistribution";

  // test for the validity of the parameter
  if ( ( p >= 0.0 ) && ( p < 1.0 ) ) {
    _p = p;
    // pre-compute the mean
    _mean = _p / (1.0 - _p);
  }
  else if ( p == 1.0 ) {
    _p = p;
    // pre-compute the mean
    _mean = INFINITE_DURATION;
  }
  else {
    // there should be an exception thrown there
    fprintf( stderr, "Error in geometricDistribution: parameter value %f",
	     p );
    fprintf( stderr, " out of bounds [0,1)\n" );
  }
}

/**
 * Calculation of the mean. Returns the value since it is pre-computed
 * @author	Alain Jean-Marie
 * @return	the mathematical expectation of the distribution
 */
double geometricDistribution::mean() {

  return _mean;

}

/**
 * Calculation of the rate, which is the inverse of the mean
 * @author	Alain Jean-Marie
 * @return	the rate
 */
double geometricDistribution::rate() {

  double res;

  if ( _mean > 0 ) {
    res = 1.0 / _mean;
  }
  else {
    res = INFINITE_RATE;
  }

  return res;

}

/**
 * Calculation of the of order n
 * @author	Alain Jean-Marie
 * @param order	order of the moment
 * @return	the moment
 */
double geometricDistribution::moment( int order ) {

  double res = 0.0;

  if ( _p == 1.0 ) {
    res = INFINITE_DURATION;
  }
  else {
    switch( order ) {
    case 0:
      res = 1.0;
      break;

    case 1:
      res = _mean;
      break;

    case 2:
      res = _p * ( 1.0 + _p ) / ( 1.0 - _p ) / ( 1.0 - _p );
      break;

    case 3:
      res = _p * ( 1.0 + 4.0*_p + _p*_p ) / pow( 1.0 - _p, 3.0 );
      break;

    case 4:
      res = _p * ( 1.0 + _p ) * ( 1.0 + 10.0*_p + _p*_p ) 
	/ pow( 1.0 - _p, 4.0 );
      break;

    default:
      fprintf( stderr, "Wwarning in geometricDistribution: cannot compute" );
      fprintf( stderr, "moment of order > 4 for geometric law. 0.0 assumed.\n" );
      res = 0.0;
    }
  }

  return res;
}

/**
 * Test of existence of a moment. These distributions always have one.
 * @author	Alain Jean-Marie
 * @param order	order of the moment
 * @return	true iff p < 1.0
 */
bool geometricDistribution::hasMoment( int order )
{
  return ( _p < 1.0 );
}

/**
 * Computation of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	LT(s)
 */
double geometricDistribution::laplace( double s )
{
  double res = 0.0;

  res = ( 1.0 - _p ) / ( 1.0 - _p * exp(-s) );

  return res;

}

/**
 * Computation of the derivative of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	d LT(s)/ds
 */
double geometricDistribution::dLaplace( double s )
{
  double res = 0.0;

  res = - _p * exp(-s) * ( 1.0 - _p ) 
    / ( 1.0 - _p * exp(-s) ) / ( 1.0 - _p * exp(-s) );

  return res;

}

/**
 * Computation of the cumulative density function at some real point x.
 * Based on the formula P(X >= k) = p^k for integer k
 * so that P( X <= x ) = P( X <= floor(x) ) = 1 - P( X >= floor(x) + 1 ) = 1 - p^( floor(x) + 1 )
 * @author	Alain Jean-Marie
 * @param x	value at which the CDF is computed
 * @return    	CDF(x)
 */
double geometricDistribution::cdf( double x )
{
  double res = 0.0;
  double k = floor(x);

  res = 1.0 - pow( _p, k );
  
  return res;

}

/**
 * Computation of the probability of some value
 * @author	Alain Jean-Marie
 * @param k	value at which the probability is computed
 * @return    	P( X = k )
 */
double geometricDistribution::getProba( double k )
{
  double res = 0.0;

  if ( ( k >= 0.0 ) && ( fabs( k - rint(k) ) < VALUE_TOLERANCE ) ) {
    res = ( 1.0 - _p ) * pow( _p, k );
  }
  
  return res;

}

// Printing a representation of the law
void geometricDistribution::write( FILE* out, int mode )
{

  switch ( mode ) {
  case PNED_PRINT_MODE:
    // not quite the "pned" mode, but an explicit description of the distrib
    fprintf( out, "Geometric on N with proba %8.4f:", _p );
    fprintf( out, " P(k) = %8.4f x (%8.4f)^k, k=0..+oo\n", 1.0 - _p, _p );
    break;

  default:
    fprintf( out, "G %8.4f ", _p );
      break;

  }

}

// Printing a representation of the law into a string
std::string geometricDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Geometric(0..infty," << _p << ")";

  return tmp.str();

}

// Rescaling the law X by some real factor f.
// Geometric distributions on N cannot be rescaled: just copied.
geometricDistribution* geometricDistribution::rescale( double factor )
{

  geometricDistribution* res = new geometricDistribution( _p );

  if ( factor != 1.0 ) {
    fprintf( stderr, "Warning in geometricDistribution: cannot be rescaled " );
    fprintf( stderr, "by factor %f <> 1.0. Just copied.\n", factor );
  }

  return res;

}

// Copying the law
geometricDistribution* geometricDistribution::copy()
{

  return this->rescale( 1.0 );

}

// Sampling from the law
double geometricDistribution::sample()
{
  double res = 0.0;

  if ( _p < 1.0 ) {
    res = floor( Distribution::exponential( -1.0/log(_p)) );
  }
  else if ( _p == 0.0 ) {
    res = 0.0;
  }
  else if ( _p == 1.0 ) {
    res = INFINITE_DURATION;
  }
  // assertion: no other case if creation enforces p in [0,1]

  return res;

}

