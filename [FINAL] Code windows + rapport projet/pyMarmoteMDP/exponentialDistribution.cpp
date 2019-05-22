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

#include "exponentialDistribution.h"
#include <stdlib.h>
#include <sstream>

// Constructor for a Exponential distribution.
exponentialDistribution::exponentialDistribution(double val)
{
  // the distribution is specified by its mean
  _mean = val;
  _name = "exponentialDistribution";

  // Calculation of the rate. It is done once, on the suspicion that
  // the rate of an exponential distribution is often needed.
  if ( _mean > 0 ) {
    _rate = 1.0 / _mean;
  }
  else if ( _mean == INFINITE_DURATION ) {
    _rate = 0.0;
  }
  else {
    _rate = INFINITE_RATE;
  }

}

// Calculation of the of order n
double exponentialDistribution::moment( int order ) {

  double res;

  switch ( order ) {
  case 0:
    res = 1.0;
    break;

  case 1:
    res = _mean;
    break;

  case 2:
    res = 2.0 * _mean * _mean;
    break;

  default:
    res = 1.0;
    for ( int i = 1; i <= order; i++ ) {
      res = res * _mean * (double)i;
    break;
    }

  }

  return res;

}

// Test of existence of a moment. These distributions always have one.
bool exponentialDistribution::hasMoment( int order ) 
{
  return true;
}

// Laplace transform at real points
double exponentialDistribution::laplace( double s )
{
  double res;

  if ( _mean != INFINITE_DURATION ) {
    res = 1.0 / ( 1.0 + s * _mean );
  }
  else {
    res = 0.0;
  }

  return res;

}

// derivative of the Laplace transform
double exponentialDistribution::dLaplace( double s )
{
  double res;

  if ( _mean != INFINITE_DURATION ) {
    res = - _mean / ( 1.0 + s * _mean ) / ( 1.0 + s * _mean );
  }
  else {
    res = 0.0;
  }

  return res;

}

// Computation of the cumulative density function at some real point x
double exponentialDistribution::cdf( double x )
{
  double res;

  if ( x < 0.0 ) {
    res = 0.0;
  }
  else {
    if ( _mean == 0.0 ) {
      res = 1.0;
    }
    else if ( _mean == INFINITE_DURATION ) {
      res = 1.0;
    }
    else {
      res = 1.0 - exp( -x / _mean );
    }
  }

  return res;

}

// output procedure
void exponentialDistribution::write( FILE* out, int mode )
{
  switch ( mode ) {
  default:
    fprintf( out, "Exponential distribution with mean %8.4f ", 
	     _mean );
  }

}

// converting the law into a descriptive string
std::string exponentialDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Exponential(m=" << _mean << ")";

  return tmp.str();

}

// Rescaling the law X by some real factor f
exponentialDistribution* exponentialDistribution::rescale( double factor )
{

  exponentialDistribution* res = new exponentialDistribution( _mean * factor );

  return res;

}

// copying the law
exponentialDistribution* exponentialDistribution::copy()
{

  return this->rescale( 1.0 );

}

// sampling
double exponentialDistribution::sample()
{

  return Distribution::exponential(_mean);

}
