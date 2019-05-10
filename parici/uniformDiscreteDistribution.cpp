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

#include "uniformDiscreteDistribution.h"
#include <stdlib.h>
#include <sstream>

/**
 * Constructor for a discrete uniform distribution on [a..b]
 * The mean is calculated at creation.
 * @author	Alain Jean-Marie
 * @param	valInf the value of a
 * @param	valSup the value of b
 * @return	an object of type uniformDiscreteDistribution
 */
uniformDiscreteDistribution::uniformDiscreteDistribution( int valInf, int valSup )
  : discreteDistribution( 0 )
{
  // Should check whether sup >= inf!
  _valInf = valInf;
  _valSup = valSup;
  _span = _valSup - _valInf;

  _name = "uniformDiscreteDistribution";

  // pre-compute the mean
  _mean = ( _valInf + _valSup ) / 2.0;

  // special case
  _isConstant = ( _valInf == _valSup );

}

double uniformDiscreteDistribution::getProba(double value) {

  double res = 0.0;

  if ( ( value > _valInf ) && ( value < _valSup ) ) {
    // determine if value is close to an integer
    int rval = (int)rint( value );
    if ( fabs( value - rval ) < VALUE_TOLERANCE ) {
      res = 1.0 / (1+_span);
    }
  }

  return res;
}

/**
 * @details Returns the value since the mean is pre-computed at creation time.
 */
double uniformDiscreteDistribution::mean() {

  return _mean;

}

/**
 * Calculation of the rate, which is the inverse of the mean
 */
double uniformDiscreteDistribution::rate() {

  double res;

  if ( _mean > 0 ) {
    res = 1.0 / _mean;
  }
  else {
    res = INFINITE_RATE;
  }

  return res;

}

// Calculation of the moment of order n
double uniformDiscreteDistribution::moment( int order ) {

  double res = 0.0;

  if ( _isConstant ) {
    res = pow( _valSup, order );
  }
  else {
    switch ( order ) {
    case 1:
      res = _mean;
      break;

    case 2:
      res = ( _valSup * ( _valSup+ 1 ) - ( _valInf-1 ) * _valInf ) / 2.0;
      break;

    case 3:
      res = ( _valSup * ( _valSup+ 1 ) * ( 2*_valSup+1 ) 
	      - ( _valInf-1 ) * _valInf * ( 2*_valInf-1 ) ) / 6.0;
      break;

    default:
      for ( int i = _valInf; i <= _valSup; i++ ) {
        res += pow( i, order );
      }
      res /= (_span+1);
    }
  }

  return res;

}

/**
 * @details These distributions have moments of all orders.
 * @return	true
 */
bool uniformDiscreteDistribution::hasMoment(int order )
{
  return true;
}

/**
 * Computation of the Laplace transform at some real point s
 */
double uniformDiscreteDistribution::laplace( double s )
{
  double res = 0.0;

  if ( _isConstant ) {
    res = exp( -s * _valInf );
  }
  else {
    res = ( s == 0.0 ? 1.0 :
	    exp( -s * _valInf ) *
	    ( 1.0 - exp( -s * (_span+1) ) ) / ( 1.0 - exp( -s ) ) );
  }

  return res;

}

/**
 * Computation of the derivative of the Laplace transform at some real point s
 */
double uniformDiscreteDistribution::dLaplace( double s )
{
  double res = 0.0;

  if ( _isConstant ) {
    res = - _valInf * exp( -s * _valInf );
  }
  else {
    res = ( s == 0.0 ?
	    _mean :
	    0.0 );
    if ( s != 0.0 ) {
      fprintf( stderr, "Warning in uniformDiscreteDistribution: dLaplace not implemented. 0 returned." );
    }
  }

  return res;

}

/**
 * Computation of the cumulative density function at some real point x
 */
double uniformDiscreteDistribution::cdf( double x )
{
  double res = 0.0;

  if ( x <= _valInf ) {
    res = 0.0;
  }
  else {
    res = ( x > _valSup ? 1.0 : rint( x - _valInf ) / (double)(1+_span) );
  }

  return res;

}

/**
 * Printing a representation of the law
 */
void uniformDiscreteDistribution::write( FILE* out, int mode )
{

  switch ( mode ) {
  default:
    fprintf( out, "uniform distribution on [%d..%d]", _valInf, _valSup );
  }

}

/**
 * Printing a representation of the law into a string
 */
std::string uniformDiscreteDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Uniform(" << _valInf << ".." << _valSup << ")";

  return tmp.str();

}

/**
 * Rescaling the law X by some real factor f
 */
uniformDiscreteDistribution* uniformDiscreteDistribution::rescale( double factor )
{

  uniformDiscreteDistribution* res =
    new uniformDiscreteDistribution( _valInf, _valSup );

  if ( factor != 1.0 ) {
    fprintf( stderr, "Warning in uniformDiscreteDistribution: cannot be rescaled " );
    fprintf( stderr, "by factor %f <> 1.0. Just copied.\n", factor );
  }

  return res;

}

/**
 * Copying the law 
 */
uniformDiscreteDistribution* uniformDiscreteDistribution::copy()
{

  return this->rescale( 1.0 );

}

/**
 * Sampling from the law 
 */
double uniformDiscreteDistribution::sample()
{

  return _valInf + (int)rint( _span * u_0_1() );

}

/**
 * Iid samples from the law in a table
 */
void uniformDiscreteDistribution::iidSample( int n, double* s )
{
  for ( int i = 0; i < n; i++ ) {
    s[i] = sample();
  }
}
