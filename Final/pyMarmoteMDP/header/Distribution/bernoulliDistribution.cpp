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

#include "bernoulliDistribution.h"
#include <stdlib.h>
#include <sstream>

// Constructor for a Bernoulli distribution.
bernoulliDistribution::bernoulliDistribution( double val )
  : discreteDistribution( 2 )
{
  _proba = val;
  _name = "bernoulliDistribution";

  // pre-compute the mean
  _mean = _proba;

  // create the image as a discrete distribution
  _nbVals = 2;
  // _values = (double*)malloc( 2*sizeof(double) );
  // _probas = (double*)malloc( 2*sizeof(double) );
  _values[0] = 0;
  _values[1] = 1;
  _probas[0] = 1.0 - _proba;
  _probas[1] = _proba;

}

/**
 * @brief Calculation of the mean. Returns the value since it is pre-computed
 * @author	Alain Jean-Marie
 * @return	the mathematical expectation of the distribution
 */
double bernoulliDistribution::mean() {

  return _mean;

}

/**
 * @brief Calculation of the rate, which is the inverse of the mean
 * @author	Alain Jean-Marie
 * @return	the rate
 */
double bernoulliDistribution::rate() {

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
 * @brief Calculation of the of order n
 * @author	Alain Jean-Marie
 * @param order	order of the moment
 * @return	the moment
 */
double bernoulliDistribution::moment( int order ) {

  double res = _proba;

  return res;

}

/**
 * @brief Test of existence of a moment. These distributions always have one.
 * @author	Alain Jean-Marie
 * @param order	order of the moment
 * @return	true
 */
bool bernoulliDistribution::hasMoment( int order )
{
  return true;
}

/**
 * @brief Computation of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	LT(s)
 */
double bernoulliDistribution::laplace( double s )
{
  double res = 1.0 - ( 1.0 - exp( -s ) ) * _proba;

  return res;

}

/**
 * @brief Computation of the derivative of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	d LT(s)/ds
 */
double bernoulliDistribution::dLaplace( double s )
{
  double res = - _proba * exp( -s );

  return res;

}

/**
 * @brief Computation of the cumulative density function at some real point x
 * @author	Alain Jean-Marie
 * @param x	value at which the CDF is computed
 * @return    	CDF(x)
 */
double bernoulliDistribution::cdf( double x )
{
  double res;

  if ( x > 1.0 ) {
    res = 1.0;
  }
  else if ( x <= 0.0 ) {
    res = 0.0;
  }
  else {
    res = 1.0 - _proba;
  }

  return res;

}

void bernoulliDistribution::write( FILE* out, int mode )
{
  switch ( mode ) {
    fprintf( out, "Bernoulli distribution with proba %8.4f\n", _proba );
  }

}

/**
 * @brief Printing a representation of the law into a string
 * @author	Alain Jean-Marie
 * @return a string
 */
std::string bernoulliDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Bernoulli(" << _proba << ")";

  return tmp.str();

}

// Rescaling the law X by some real factor f
bernoulliDistribution* bernoulliDistribution::rescale( double factor )
{

  bernoulliDistribution* res = new bernoulliDistribution( _proba );

  if ( factor != 1.0 ) {
    fprintf( stderr, "Warning in bernoulliDistribution: cannot be rescaled " );
    fprintf( stderr, "by factor %f <> 1.0. Just copied.\n", factor );
  }

  return res;

}

// Copying the law
bernoulliDistribution* bernoulliDistribution::copy()
{

  return this->rescale( 1.0 );

}

// Sampling from the law
double bernoulliDistribution::sample()
{
  double res;

  if ( Distribution::u_0_1() < _proba ) {
    res = 1.0;
  }
  else {
    res = 0.0;
  }

  return res;

}
