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

#include "header/discreteDistribution.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>

// creation just from the size; usable only from siblings
discreteDistribution::discreteDistribution( int sz )
{

  _name = "discreteDistribution";
  _nbVals = sz;
  _values = (double*) malloc( _nbVals * sizeof(double) );
  _probas = (double*) malloc( _nbVals * sizeof(double) );

  // pre-compute the mean, with a default 0 value
  _mean = 0.0;

}

// standard constructor from arrays
discreteDistribution::discreteDistribution( int sz, 
                                            double *vals, double *probas )
{

  _name = "discreteDistribution";
  _nbVals = sz;
  _values = (double*) malloc( _nbVals * sizeof(double) );
  _probas = (double*) malloc( _nbVals * sizeof(double) );
  for ( int i = 0; i < _nbVals; i++ ) {
    _values[i] = vals[i];
    _probas[i] = probas[i];
  }

  // pre-compute the mean
  _mean = 0.0;
  for ( int i=0; i < _nbVals; i++ ) {
    _mean += _values[i] * _probas[i];
  }

}

// constructor from a file
discreteDistribution::discreteDistribution( int sz, char *name )
{

  _name = "discreteDistribution";
  _nbVals = sz;
  _values = (double*) malloc( _nbVals * sizeof(double) );
  _probas = (double*) calloc( _nbVals, sizeof(double) );
  for ( int i = 0; i < _nbVals; i++ ) {
    _values[i] = i;
  }

  FILE *in = fopen( name, "r" );
  int res;
  if ( NULL != in ) {
    for ( int i = 0; i < _nbVals; i++ ) {
      res = fscanf( in, "%lf\n", &(_probas[i]) );
    }
    if ( EOF == res ) {
      fprintf( stderr, "Warning: short file '%s'. Some zero probas assumed.\n",
	       name );
    }
  }
  else {
    fprintf( stderr, 
	     "Warning: cannot open file '%s' for reading. Zero probas assumed.\n",
	     name );
  }
  fclose( in );

  // pre-compute the mean
  _mean = 0.0;
  for ( int i=0; i < _nbVals; i++ ) {
    _mean += _values[i] * _probas[i];
  }

}

// standard destructor
discreteDistribution::~discreteDistribution()
{

  if ( (double*)NULL != _values ) free( _values );
  if ( (double*)NULL != _probas ) free( _probas );

}

// accessor for elements of the probas array
double discreteDistribution::getProbaByIndex(int i)
{
  double res = 0.0;

  if ( ( i >= 0 ) && ( i < _nbVals ) ) {
    res = _probas[i];
  }
  else {
    fprintf( stdout, "Warning: accessing entry '%d' out of range [0..%d]",
	     i, _nbVals );
    fprintf( stdout, " in discreteDistribution object. 0.0 assumed.\n" );
  }

  return res;

}

// computation of the probability distribution
double discreteDistribution::getProba(double value)
{
  double res = 0.0;

  for ( int i = 0; i < _nbVals; i++ ) {
    if ( fabs( _values[i] - value ) < VALUE_TOLERANCE ) {
      res = _probas[i];
    }
  }

  return res;

}
// accessor for values
double discreteDistribution::getValue(int i)
{
  double res = 0.0;

  if ( ( i >= 0 ) && ( i < _nbVals ) ) {
    res = _values[i];
  }
  else {
    fprintf( stdout, "Warning: accessing entry '%d' out of range [0..%d]",
	     i, _nbVals );
    fprintf( stdout, " in discreteDistribution object. 0.0 assumed.\n" );
  }

  return res;

}

// write accessor for probas
bool discreteDistribution::setProba(int i,double v)
{
  bool res = true;

  if ( ( i >= 0 ) && ( i < _nbVals ) ) {
    _probas[i] = v;
  }
  else {
    fprintf( stdout, "Warning: try to set entry '%d' out of range [0..%d]",
	     i, _nbVals );
    fprintf( stdout, " in discreteDistribution object.\n" );
    res = false;
  }

  return res;

}

// the rate is the inverse of the mean
double discreteDistribution::rate() {

  double res;

  if ( _mean > 0 ) {
    res = 1.0 / _mean;
  }
  else {
    res = INFINITE_RATE;
  }

  return res;

}

// moments
double discreteDistribution::moment( int order ) {

  double res = 0.0;

  for ( int i=0; i < _nbVals; i++ ) {
    res += _probas[i] * pow( _values[i], order );
  }

  return res;

}

// test for moments
bool discreteDistribution::hasMoment( int order ) 
{
  return true;
}

/**
 * @brief Computation of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	LT(s)
 */
double discreteDistribution::laplace( double s )
{
  double res = 0.0;

  for ( int i=0; i < _nbVals; i++ ) {
    res += exp( -s * _values[i] ) * _probas[i];
  }

  return res;

}

/**
 * @brief Computation of the derivative of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	d LT(s)/ds
 */
double discreteDistribution::dLaplace( double s )
{
  double res = 0.0;

  for ( int i=0; i < _nbVals; i++ ) {
    res -= _values[i] * exp( -s * _values[i] ) * _probas[i];
  }

  return res;

}

/**
 * @brief Computation of the cumulative density function at some real point x
 * @author	Alain Jean-Marie
 * @param x	value at which the CDF is computed
 * @return    	CDF(x)
 */
double discreteDistribution::cdf( double x )
{
  double res = 0.0;

  for ( int i=0; i < _nbVals; i++ ) {
    if ( x > _values[i] )
      res += _probas[i];
  }

  return res;

}

/**
 * @brief Printing a representation of the law
 * @author	Alain Jean-Marie
 * @param out file descriptor of the output stream
 * @param mode	representation for the output
 */
void discreteDistribution::write( FILE* out, int mode )
{
  switch ( mode ) {
  case PNED_PRINT_MODE:
    fprintf( out, "discrete values { ");
    for ( int i=0; i < _nbVals; i++ )
      fprintf( out, "%8.4f ", _values[i] );
    fprintf( out, "} probas { ");
    for ( int i=0; i < _nbVals; i++ )
      fprintf( out, "%8.4f ", _probas[i] );
    fprintf( out, "} ");
    break;

  case MAPLE_PRINT_MODE:
    fprintf( out, "Vector( [" );
    for ( int i=0; i < _nbVals-1 ; i++ )
      fprintf( out, "  %12e,\n ", _probas[i] );
    fprintf( out, "  %12e\n", _probas[_nbVals-1] );
    fprintf( out, "] );\n" );
    break;

  default:
    fprintf( out, "discrete [ ");
    for ( int i=0; i < _nbVals; i++ )
      fprintf( out, "%8.4f ", _values[i] );
    fprintf( out, "] [ ");
    for ( int i=0; i < _nbVals; i++ )
      fprintf( out, "%8.4f ", _probas[i] );
    fprintf( out, "] ");
  }

}

/**
 * @brief Printing a representation of the law into a string
 * @author	Alain Jean-Marie
 * @return a string
 */
std::string discreteDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Discrete on (";
  for ( int i=0; i < _nbVals-1; i++ ) {
    tmp << _values[i] << ",";
  }
  tmp << _values[_nbVals-1] << ") wp (";
  for ( int i=0; i < _nbVals-1; i++ ) {
    tmp << _probas[i] << ",";
  }
  tmp << _probas[_nbVals-1] << ")";

  return tmp.str();

}

/**
 * @brief Rescaling the law X by some real factor f
 * @author	Alain Jean-Marie
 * @param factor the factor f
 * @result	the law f*X
 */
discreteDistribution* discreteDistribution::rescale( double factor )
{
  double *newVals = (double*) malloc( _nbVals * sizeof(double) );

  for ( int i=0; i < _nbVals; i++ ) {
    newVals[i] = factor * _values[i];
  }

  discreteDistribution* res = new discreteDistribution( _nbVals, newVals, _probas );

  free( newVals );

  return res;

}

/**
 * @brief Copying the law
 * @author	Alain Jean-Marie
 * @result	a copy of the law
 */
discreteDistribution* discreteDistribution::copy()
{

  return this->rescale( 1.0 );

}

/**
 * @brief Sampling from the law
 * This is the straightforward, non optimized, linear-time algorithm
 * @author	Alain Jean-Marie
 * @result	a copy of the law
 */
double discreteDistribution::sample()
{
  register double	u;
  int			i;

  u = u_0_1();

  for ( i=0; (i<_nbVals) && (u>_probas[i]); i++ )
    u -= _probas[i];

  return _values[i];
}

/**
 * @brief Computes the L1 distance between this distribution and another one.
 * In case of incompatible or infinite sizes, a negative number is returned.
 * @author	Alain Jean-Marie
 * @param d	the second distribution
 * @return 	the distance
 */
double discreteDistribution::distanceL1(discreteDistribution* d)
{
  double res;

  // assertion: all discreteDistribution objects has a **finite** state space
  int m1 = _nbVals;
  int m2 = d->nbVals();
  // arranging so that d1 has the longest list of values
  int mmin = ( m1 < m2 ? m1 : m2 );
  int mmax = ( m1 > m2 ? m1 : m2 );
  discreteDistribution *d1;
  discreteDistribution *d2;
  if ( m2 > m1 ) {
    d1 = d;
    d2 = this;
  }
  else {
    d1 = this;
    d2 = d;
  }

  res = 0.0;
  for ( int i=0; i < mmin; i++ ) {
    res += fabs( d1->getProba(i) - d2->getProba(i) );
  }
  for ( int i=mmin; i < mmax; i++ ) {
    res += d1->getProba(i);
  }

  // else {
  //   std::cerr << "Error in distanceL1: "
  // 	      << "incompatible state space size (" << d->nbVals()
  // 	      << "<>" << _nbVals << ")" << std::endl;
  //   //The method should throw an exception/error
  //   res = -1.0; // a negative number signals an error!
  // }

  return res;
}

/**
 * @brief Computes the L2 distance between two distributions.
 * In case of incompatible or infinite sizes, a negative number is returned.
 * @author	Alain Jean-Marie
 * @param d	the second distribution
 * @return 	the distance
 */
double discreteDistribution::distanceL2( discreteDistribution* d )
{
  double res;

  if ( ( d->nbVals() == _nbVals ) || ( _nbVals != INFINITE_STATE_SPACE_SIZE ) ) {
    res = 0.0;
    for ( int i=0; i < _nbVals; i++ ) {
      res += ( d->getProba(i) - _probas[i] )*( d->getProba(i) - _probas[i] );
    }
    res = sqrt( res );
  }
  else {
    fprintf( stderr, "Error in Distrib_Distance_L2: " );
    fprintf( stderr, "incompatible state space size (%d<>%d)\n",
         d->nbVals(), _nbVals );
	/* The method should throw an exception/error */
    res = -1.0; /* a negative number signals an error! */
  }
  
  return res;
}

/**
 * @brief Computes the L-infinity distance between two distributions.
 * In case of incompatible or infinite sizes, a negative number is returned.
 * @author	Alain Jean-Marie
 * @param d	the second distribution
 * @return 	the distance
 */
double discreteDistribution::distanceLinfinity( discreteDistribution* d )
{
  double res;
  double val;

  if ( ( d->nbVals() == _nbVals ) || ( _nbVals != INFINITE_STATE_SPACE_SIZE ) ) {
    res = 0.0;
    for ( int i=0; i < _nbVals; i++ ) {
      res = ( val = fabs( d->getProba(i) - _probas[i] ) > res ? val : res );
    }
  }
  else {
    fprintf( stderr, "Error in distanceLinfinity: " );
    fprintf( stderr, "incompatible state space size (%d<>%d)\n",
         d->nbVals(), _nbVals );
	/* The method should throw an exception/error */
    res = -1.0; /* a negative number signals an error! */
  }
  
  return res;
}
