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

#include "uniformDistribution.h"
#include <stdlib.h>
#include <sstream>

/**
 * Constructor for a continuous uniform distribution on [a,b]
 * The mean is calculated at creation.
 * @author	Alain Jean-Marie
 * @param	inf the value of a
 * @param	sup the value of b
 * @return	an object of type uniformDistribution
 */
uniformDistribution::uniformDistribution( double inf, double sup )
{
  // Should check whether sup >= inf!
  _valInf = inf;
  _valSup = sup;
  _span = _valSup - _valInf;

  _name = "uniformDistribution";

  // pre-compute the mean
  _mean = ( _valInf + _valSup ) / 2.0;

  // special case
  _isConstant = ( _valInf == _valSup );

}

/**
 * @details Returns the value since the mean is pre-computed at creation time.
 */
double uniformDistribution::mean() {

  return _mean;

}

/**
 * Calculation of the rate, which is the inverse of the mean
 */
double uniformDistribution::rate() {

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
 */
double uniformDistribution::moment(int order )
{
  double res;

  if ( _isConstant ) {
    res = pow( _valSup, order );
  }
  else {
    res = pow( _valSup, 1+order ) - pow( _valInf, 1+order );
    res /= _span;
    res /= ( 1.0 + (double)order );
  }

  return res;

}

/**
 * @details These distributions have moments of all orders.
 * @return	true
 */
bool uniformDistribution::hasMoment(int order )
{
  return true;
}

/**
 * Computation of the Laplace transform at some real point s
 */
double uniformDistribution::laplace( double s )
{
  double res;

  if ( _isConstant ) {
    res = exp( -s * _valInf );
  }
  else {
    res = ( s == 0.0 ? 1.0 :
	    ( exp( -s * _valInf ) - exp( -s * _valSup ) )
	    / s / _span );
  }

  return res;

}

/**
 * Computation of the derivative of the Laplace transform at some real point s
 */
double uniformDistribution::dLaplace( double s )
{
  double res = 0.0;

  if ( _isConstant ) {
    res = - _valInf * exp( -s * _valInf );
  }
  else {
    res = ( s == 0.0 ?
	    _span / 2.0 :
	    ( 1.0 + s * _valInf ) * exp( -s * _valInf ) - 
	    ( 1.0 + s * _valSup ) * exp( -s * _valSup ) / s / s /
	    _span );
  }

  return res;

}

/**
 * Computation of the cumulative density function at some real point x
 */
double uniformDistribution::cdf( double x )
{
  double res = 0.0;

  if ( x <= _valInf ) {
    res = 0.0;
  }
  else {
    res = ( x > _valSup ? 1.0 : ( x - _valInf ) / _span );
  }

  return res;

}

/**
 * Printing a representation of the law
 */
void uniformDistribution::write( FILE* out, int mode )
{

  switch ( mode ) {
  default:
    fprintf( out, "uniform distribution on [%8.4f,%8.4f]", _valInf, _valSup );
  }

}

/**
 * Printing a representation of the law into a string
 */
std::string uniformDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Uniform(" << _valInf << "," << _valSup << ")";

  return tmp.str();

}

/**
 * Rescaling the law X by some real factor f
 */
uniformDistribution* uniformDistribution::rescale( double factor )
{

  uniformDistribution* res = new uniformDistribution( factor*_valInf,
                                                      factor*_valSup );

  return res;

}

/**
 * Copying the law 
 */
uniformDistribution* uniformDistribution::copy()
{

  return this->rescale( 1.0 );

}

/**
 * Sampling from the law 
 */
double uniformDistribution::sample()
{

  return _valInf + _span * u_0_1();

}

/**
 * Iid samples from the law in a table
 */
void uniformDistribution::iidSample( int n, double* s )
{
  for ( int i = 0; i < n; i++ ) {
    s[i] = sample();
  }
}
