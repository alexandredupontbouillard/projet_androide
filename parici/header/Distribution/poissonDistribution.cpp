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

#include "poissonDistribution.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>

#ifdef WITH_R
RInside* poissonDistribution::_Rmotor = (RInside*)NULL;
#else
typedef void* RInside;
#endif

poissonDistribution::poissonDistribution( double lambda )
{

  _name = "poissonDistribution";
  _lambda = lambda;

  // pre-compute the mean
  _mean = lambda;

}

poissonDistribution::~poissonDistribution()
{
  // Useless attempt at recycling the RInside engine.
  // Once RInside is initialized, it stays so even when the object
  // is destroyed.
// #ifdef WITH_R
//  delete _Rmotor;
//  _Rmotor = (RInside*)NULL;
// #endif
}

RInside* poissonDistribution::Rmotor()
{
#ifdef WITH_R
  if ( _Rmotor == (RInside*)NULL ) {
    _Rmotor = new RInside(0,NULL);
    _Rmotor->setVerbose(true);
  }
  return _Rmotor;
#else
  try{
    throw "Error: functionality available only when compiled with R";
  }
  catch (char const* e)
  {
    std::cout << "An exception occurred.  " << e <<'\n';
  }
#endif
  return NULL;
}

/**
 * Calculation of the mean. Returns the value since it is pre-computed
 * @author	Alain Jean-Marie
 * @return	the mathematical expectation of the distribution
 */
double poissonDistribution::mean() {

  return _mean;

}

/**
 * Calculation of the rate, which is the inverse of the mean
 * @author	Alain Jean-Marie
 * @return	the rate
 */
double poissonDistribution::rate() {

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
double poissonDistribution::moment( int order ) {

  double res = _mean; // default value which is moments of order 1 and 2

  switch( order ) {
  case 0:
    res = 1.0;
    break;

  case 1:
    res = _mean;
    break;

  case 2:
    res = _lambda * ( _lambda + 1 );
    break;

  case 3:
    res = _lambda * ( _lambda*_lambda + 3*_lambda + 1 );
    break;

  case 4:
    res = _lambda * ( _lambda*_lambda*_lambda + 6*_lambda*_lambda + 7*_lambda + 1 );
    break;

  default:
    fprintf( stderr, "Warning in poissonDistribution: cannot compute" );
    fprintf( stderr, "moment of order > 4 for Poisson law. 0.0 assumed.\n" );
    res = 0.0;
  }

  return res;

}

// Test of existence of a moment. These distributions always have one... or not.
bool poissonDistribution::hasMoment( int order )
{
  return true;
}

/**
 * Computation of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	LT(s)
 */
double poissonDistribution::laplace( double s )
{
  double res = 0.0;

  res = exp( _lambda * ( exp(-s) - 1 ) );

  return res;

}

/**
 * Computation of the derivative of the Laplace transform at some real point s
 * @author	Alain Jean-Marie
 * @param s	value at which the derivative is computed
 * @return	d LT(s)/ds
 */
double poissonDistribution::dLaplace( double s )
{
  double res = 0.0;

  res = - _lambda * exp(-s) * exp( _lambda * ( exp(-s) - 1 ) );

  return res;

}

/**
 * Computation of the cumulative density function at some real point x
 * @author	Alain Jean-Marie
 * @param x	value at which the CDF is computed
 * @return    	CDF(x)
 */
double poissonDistribution::cdf( double x )
{
  double res = 0.0;
  double incr = exp( -_lambda );

  // straightforward method, probably not ideal for large x and/or large lambda
  for ( int i = 0; i <= x; i++ ) {
    res += incr;
    incr *= _lambda/(i+1);
  }

  return res;

}

/**
 * Computation of the probability of some value
 * @author	Alain Jean-Marie
 * @param k	value at which the probability is computed
 * @return    	P( X = k )
 */
double poissonDistribution::getProba( double k )
{
  double res = 0.0;

  if ( ( k >= 0.0 ) && ( fabs( k - rint(k) ) < VALUE_TOLERANCE ) ) {
    // use logarithmic method to improve stability
    double logres = -_lambda;
    double loglam = log( _lambda );

    for ( int i = 1; i <= k; i++ ) {
      logres += loglam - log(i);
    }

    res = exp(logres);
  }

  return res;

}

/**
 * Printing a representation of the law
 * @author	Alain Jean-Marie
 * @param out the output stream
 * @param mode	representation for the output
 */
void poissonDistribution::write( FILE* out, int mode )
{

  switch ( mode ) {
  case PNED_PRINT_MODE:
    // not quite the "pned" mode, but an explicit description of the distrib
    fprintf( out, "Poisson with rate %8.4f:", _lambda );
    fprintf( out, " P(k) = exp(-%8.4f) (%8.4f)^k / k!, k=0..+oo\n", _lambda, _lambda );
    break;

  default:
    fprintf( out, "P %8.4f ", _lambda );
      break;

  }

}

/**
 * Printing a representation of the law into a string
 * @author	Alain Jean-Marie
 * @return a string
 */
std::string poissonDistribution::toString()
{
  std::string res;
  std::ostringstream tmp;

  tmp << "Poisson(" << _lambda << ")";

  return tmp.str();

}

/**
 * Rescaling the law X by some real factor f
 * @author	Alain Jean-Marie
 * @param factor the factor f
 * @result	the law f*X
 */
poissonDistribution* poissonDistribution::rescale( double factor )
{

  poissonDistribution* res = new poissonDistribution( _lambda * factor );

  return res;

}

/**
 * Copying the law 
 * @author	Alain Jean-Marie
 * @result	a copy of the law
 */
poissonDistribution* poissonDistribution::copy()
{

  return this->rescale( 1.0 );

}

// Sampling from the law. Call to the "R" function, if available.
double poissonDistribution::sample()
{
  // register double	u;
  double res = 0.0;

#ifdef WITH_R
  // retreive R engine, with creation if needed.
  _Rmotor = Rmotor();
  std::ostringstream tmp;
  tmp << "as.matrix( rpois(1," << _lambda << ") )";
  std::string cmd = tmp.str();
  // It appears that the matrix has only one row and is handled as a vector
  Rcpp::NumericMatrix M = ((SEXP) _Rmotor->parseEval(cmd));
  res = M[0];
  // const char* format = (Rcpp::as<std::string>(val)).c_str();
  // if ( 1 != sscanf( format, "%lf", &res ) ) {
  //  fprintf( stderr,
  //           "Internal error in %s: '%s' does not contain a number. 0 assumed.\n",
  //           "poissonDistribution::sample()", format );
  //}
  // to be finished
#else
  fprintf( stderr, "Warning in poissonDistribution::sample(): only available with R. 0.0 returned.\n");
#endif

  return res;

}

// iid samples
void poissonDistribution::iidSample( int n, double* s )
{
#ifdef WITH_R
  // retreive R engine, with creation if needed.
  _Rmotor = Rmotor();
  std::ostringstream tmp ;
  tmp << "as.matrix( rpois(" << n << "," << _lambda << ") )";
  std::string cmd = tmp.str();
  Rcpp::NumericMatrix M = ((SEXP) _Rmotor->parseEval(cmd));
  // The result appears to be a column vector
  if ( n != M.nrow() ) {
    fprintf( stderr, "Warning in poissonDistribution::iidSample(): %d samples got, %d expected.\n",
             M.nrow(), n );
  }
  for ( int i = 0; ( i < n ) && ( i < M.nrow() ); i++ ) {
    s[i] = M[i];
  }
#else
  for ( int i = 0; i < n; i++ ) {
    s[i] = sample();
  }
#endif
}
