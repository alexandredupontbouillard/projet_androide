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

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>
#include "marmoteConstants.h"

// Definition of +oo for durations and rates
#define INFINITE_DURATION       HUGE_VAL
#define INFINITE_RATE	        HUGE_VAL

// Definitions for printing styles
#define NO_PRINT_MODE		0
#define	DEFAULT_PRINT_MODE	1
#define	NORMAL_PRINT_MODE	1
#define QNAP_PRINT_MODE		2
#define PROSIT_PRINT_MODE	3
#define MMPP_PRINT_MODE		4
#define PNED_PRINT_MODE		5
#define MAPLE_PRINT_MODE    6

// Definition of the class
/**
 * @brief A class for representing probability distributions.
 *
 */
class Distribution {

 protected:
  std::string _name; /**< a string uniquely representing the distribution */

 protected:
  // Distribution();
  // Distribution( Distribution d ); // copy

 public:
  /**
   * @brief Standard destructor
   */
  virtual ~Distribution() {}

public:
  // accessors
  /**
   * @brief Read accessor to the type name of the distribution
   *
   * @return the name
   */
  std::string name() { return _name; }

 public:
  // probabilistic member functions
  /**
   * @brief computing the mathematical expectation or mean
   *
   * @return the mean
   */
  virtual double	mean() = 0;
  /**
   * @brief computing the "rate", defined as the inverse of the mean
   *
   * @return the rate
   */
  virtual double	rate() = 0;
  /**
   * @brief Computing the moments of the distribution
   *
   * @param order the order of the moment to be computed
   * @return the moment
   */
  virtual double	moment( int order ) = 0;
  /**
   * @brief Computing the variance of the random variable: the second moment minus
   * the square of the first moment. Variance is the square of the coefficient of
   * variation.
   * The Distribution class offers a default implementation.
   *
   * @return the variance
   */
  double		variance();
  /**
   * @brief computing the Laplace transform of the distribution at real point
   *
   * @param s the value at which to compute
   * @return the value of the Laplace transform
   */
  virtual double	laplace( double s ) = 0;
  /**
   * @brief computing the derivative of the Laplace transform at real points
   *
   * @param s the value at which to compute
   * @return the value of the derivative of the Laplace transform
   */
  virtual double	dLaplace( double s ) = 0;
  /**
   * @brief computing the cumulative distribution function at some real point x.
   * This is the probability that the random variable is less or equal to x.
   *
   * @param x the value at which to compute the cdf
   * @return the value of the cdf
   */
  virtual double	cdf( double x ) = 0;
  /**
   * @brief computing the complementary cumulative distributon function (or tail)
   * at some real point x. This is the probability that the random variable is
   * strictly larger than x. 
   * The Distribution class offers a default implementation.
   *
   * @param x the value at which to compute the ccdf
   * @return the value of the ccdf
   */
  double	ccdf( double x ) { return 1.0 - cdf(x); };
  /**
   * @brief test for the existence of moments of any order
   *
   * @param order the order of the moment to be tested
   * @return true if the moment exists, false otherwise
   */
  virtual bool		hasMoment( int order ) = 0;
  
  /**
   * @brief rescaling a distribution by some real factor. Not all distributions
   * allow this for any real factor. If the operation fails, or if the factor
   * is 1.0, a copy of the 
   * distribution should be returned (<b>not</b> by using the copy() function).
   * @see copy()
   *
   * @param factor the real factor to be used
   * @return a new distribution object
   */
  virtual Distribution* rescale( double factor ) = 0;
  /**
   * @brief copying a distribution. Typically implemented as rescale(1.0).
   * @see rescale(double)
   *
   * @return a copy of the distribution
   */
  virtual Distribution* copy() = 0;
  /**
   * @brief drawing a (pseudo)random value according to the distribution.
   *
   * @return a sample
   */
  virtual double	sample() = 0;
  /**
   * @brief drawing an i.i.d. sample from the distribution. The result is returned
   * in an array (that must have been already allocated) passed as a parameter.
   * The Distribution class offers the default implementation with repeated call to
   * sample().
   *
   * @param n the number of values to sample
   * @param s an array to be filled with the sample
   */
  void		iidSample( int n, double* s );

  /**
   * @brief Computing generally the L1 distance between distributions.
   * @author Alain Jean-Marie
   * @param d the distribution with which the distance is computed
   * @return the distance
   */
  virtual double	distanceL1( Distribution* d );

  /**
   * @brief Property test function. Current properties are:
   * - discreteDistribution: tests if the distribution is discrete and finite
   * - discrete: tests if the distribution is discrete
   * - continuous: tests if the distribution is continuous.
   * @author Alain Jean-Marie
   * @param pro a string/key for the property to be tested
   * @return true if the distribution has the property, false otherwise
   */
  virtual bool		hasProperty( std::string pro );
  /* functions existing in the ERS package, not implemented yet
  Law_Desc Parse_Law( char Nom, Liste_Reel Params, 
		    Liste_Reel Params2, Law_List The_SubLaws );
  int Parse_Law_From_Args ( char**, int, int, Law_Desc );
  Law_List Append_Law( Law_List, Law_Desc New_Law );
  */

 public:
  // standard pseudo-random generators
  /**
   * The generator for the uniform distribution on [0,1]. Serves as the
   * building block for most sampling algorithms.
   * @author Alain Jean-Marie
   * @return a sample of the Uniform([0,1]) distribution
   */
  static double u_0_1(void);
  /**
   * The generator for exponential distributions. Warning: the parameter
   * is the mean.
   * @author Alain Jean-Marie
   * @param mean the average (not the rate)
   * @return a sample of the Exponential distribution with the given mean
   */
  static double exponential(double mean);

 public:
  /**
   * @brief an utility to convert the distribution into a string.
   *
   * @return the string representing the distribution.
   */
  virtual std::string toString() = 0;
  /**
   * @brief an utility to write the distribution to some file, according 
   * to some format.
   *
   * @param out the file descriptor of the file
   * @param mode a code for the format to be used
   */
  virtual void write( FILE *out, int mode ) = 0;
  //  virtual void fprint() = 0; 
  /**
   * @brief write on stdout with NORMAL_PRINT_MODE
   *
   */
  void fprint() { write( stdout, NORMAL_PRINT_MODE ); }

 protected:
  static const double VALUE_TOLERANCE; /**< tolerance accepted when looking for values */
  double _mean; /**< the mathematical expectation of the distribution */

};

#endif // DISTRIBUTION_H
