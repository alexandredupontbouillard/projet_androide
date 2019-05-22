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

#include "templateDistribution.h"
#include <stdlib.h>
#include <sstream>

templateDistribution::templateDistribution( int x, double y, double* z )
{

  _name = "templateDistribution";

  // pre-compute the mean
  _mean = 0.0;

}

/**
 * Calculation of the mean. Returns the value since it is pre-computed
 * @author	Alain Jean-Marie
 * @return	the mathematical expectation of the distribution
 */
double templateDistribution::mean() {

  return _mean;

}

/**
 * Calculation of the rate, which is the inverse of the mean
 * @author	Alain Jean-Marie
 * @return	the rate
 */
double templateDistribution::rate() {

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
double templateDistribution::moment(int order ) {

  double res = 0.0;

  return res;

}

// Test of existence of a moment. These distributions always have one... or not.
bool templateDistribution::hasMoment(int order )
{
  return true;
}

/**
 * Computation of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	LT(s)
 */
double templateDistribution::laplace( double s )
{
  double res = 0.0;

  return res;

}

/**
 * Computation of the derivative of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	d LT(s)/ds
 */
double templateDistribution::dLaplace( double s )
{
  double res = 0.0;

  return res;

}

/**
 * Computation of the cumulative density function at some real point x
 * @author	Alain Jean-Marie
 * @param x	value at which the CDF is computed
 * @return    	CDF(x)
 */
double templateDistribution::cdf( double x )
{
  double res = 0.0;

  return res;

}

/**
 * Printing a representation of the law
 * @author	Alain Jean-Marie
 * @param out the output stream
 * @param mode	representation for the output
 */
void templateDistribution::write( FILE* out, int mode )
{

  switch ( mode ) {
  case PNED_PRINT_MODE:

    break;

  default:

  }

}

/**
 * Printing a representation of the law into a string
 * @author	Alain Jean-Marie
 * @return a string
 */
std::string templateDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Template(" << ... << ")";

  return tmp.str();

}

/**
 * Rescaling the law X by some real factor f
 * @author	Alain Jean-Marie
 * @param factor the factor f
 * @result	the law f*X
 */
templateDistribution* templateDistribution::rescale( double factor )
{

  templateDistribution* res = new templateDistribution( ... );

  return res;

}

/**
 * Copying the law 
 * @author	Alain Jean-Marie
 * @result	a copy of the law
 */
templateDistribution* templateDistribution::copy()
{

  return this->rescale( 1.0 );

}

/**
 * Sampling from the law 
 * This is the straightforward, non optimized, linear-time algorithm
 * @author	Alain Jean-Marie
 * @result	a copy of the law
 */
double templateDistribution::sample()
{
  register double	u;

  u = u_0_1();

}

// iid samples
void templateDistribution::iidSample( int n, double* s )
{
  for ( int i = 0; i < n; i++ ) {
    s[i] = sample();
  }
}
