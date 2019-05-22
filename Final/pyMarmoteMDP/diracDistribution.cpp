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

#include "header/diracDistribution.h"
#include <stdlib.h>
#include <sstream>

/**
 * Constructor for a Dirac distribution.
 * The mean is calculated at creation.
 * @author	Alain Jean-Marie
 * @param	val the value of the deterministic law
 * @return	an object of type diracDistribution
 */
diracDistribution::diracDistribution( double val )
  : discreteDistribution( 1 )
{
  _value = val;
  _name = "diracDistribution";

  // pre-compute the mean
  _mean = _value;

  // create the image as a discrete distribution
  _nbVals = 1;
  // _values = (double*)malloc( sizeof(double) );
  // _probas = (double*)malloc( sizeof(double) );
  _values[0] = _value;
  _probas[0] = 1.0;

}

// computing the proba
double diracDistribution::getProba(double value)
{
  double res = 0.0;

  if ( fabs( value - _value) < VALUE_TOLERANCE ) {
    res = 1.0;
  }

  return res;
}

/**
 * Calculation of the mean. Returns the value since it is pre-computed
 * @author	Alain Jean-Marie
 * @return	the mathematical expectation of the distribution
 */
double diracDistribution::mean() {

  return _mean;

}


/**
 * Calculation of the rate, which is the inverse of the mean
 * @author	Alain Jean-Marie
 * @return	the rate
 */
double diracDistribution::rate() {

  double res;

  if ( _mean > 0 ) {
    res = 1.0 / _mean;
  }
  else {
    res = INFINITE_RATE;
  }

  return res;

}

// Calculation of the of order n
double diracDistribution::moment( int order ) {

  double res = pow( _value, order );

  return res;

}

// calculation of the variance
double diracDistribution::variance()
{
  return 0.0;
}

/**
 * Test of existence of a moment. These distributions always have one.
 * @author	Alain Jean-Marie
 * @param order	order of the moment
 * @return	true
 */
bool diracDistribution::hasMoment( int order )
{
  return true;
}

/**
 * Computation of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	LT(s)
 */
double diracDistribution::laplace( double s )
{
  double res = exp( -s * _value);

  return res;

}

/**
 * Computation of the derivative of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	d LT(s)/ds
 */
double diracDistribution::dLaplace( double s )
{
  double res = - _value * exp( -s * _value );

  return res;

}

/**
 * Computation of the cumulative density function at some real point x
 * @author	Alain Jean-Marie
 * @param x	value at which the CDF is computed
 * @return    	CDF(x)
 */
double diracDistribution::cdf( double x )
{
  double res;

  if ( x > _value ) {
    res = 1.0;
  }
  else {
    res = 0.0;
  }

  return res;

}

void diracDistribution::write( FILE* out, int mode )
{

  switch ( mode ) {
  default:
    fprintf( out, "Dirac distribution at %f\n", _value );
  }

}

/**
 * Printing a representation of the law into a string
 * @author	Alain Jean-Marie
 * @return a string
 */
std::string diracDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Dirac(" << _value << ")";

  return tmp.str();

}

// Rescaling the law X by some real factor f
diracDistribution* diracDistribution::rescale( double factor )
{

  diracDistribution* res = new diracDistribution( _value * factor );

  return res;

}

// Copying the law
diracDistribution* diracDistribution::copy()
{

  return this->rescale( 1.0 );

}

// Sampling from the law
double diracDistribution::sample()
{
  return _value;

}

// Iid samples from the law in a table
void diracDistribution::iidSample( int nbSamples, double* sequence )
{
  for ( int i = 0; i < nbSamples; i++ ) {
    sequence[i] = _value;
  }
}
