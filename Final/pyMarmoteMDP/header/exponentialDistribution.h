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

#ifndef EXPONENTIALDISTRIBUTION_H
#define EXPONENTIALDISTRIBUTION_H

#include "discreteDistribution.h"

/**
 * @brief The class representing the (negative) exponential distribution
 *
 */
class exponentialDistribution : public virtual Distribution {

public:
  // constructors
/**
 * @brief Unique constructor for the exponential distribution, from its average.
 * The rate is computed at creation time.
 * @author Alain Jean-Marie
 * @param val the mean of the distribution
 */
  exponentialDistribution( double val ); // creation from a value

private:
  // specific variables
  double _rate; /**< the rate is the inverse of the mean */

public:
  // accessors
  /**
   * @brief Read accessor for the rate.
   *
   * @return the value of the rate
   */
  double	rate() { return _rate; }

public:
  // probabilistic member functions
  /**
   * @brief Calculation of the mean. Returns the value since it is pre-computed
   * @author	Alain Jean-Marie
   * @return	the mathematical expectation of the distribution
   */
  double	mean() { return _mean; };

  /**
   * @copydoc Distribution::moment()
   */
  double	moment( int order );
  /**
   * @copydoc Distribution::laplace(double s)
   */
  double	laplace( double s ); // Laplace transform at real points
  /**
   * @copydoc Distribution::dLaplace(double s)
   */
  double	dLaplace( double s ); // derivative of the Laplace transform
  /**
   * @copydoc Distribution::cdf(double x)
   */
  double	cdf( double x );
  /**
   * @brief Test of existence of a moment. These distributions always have one.
   * @author	Alain Jean-Marie
   * @param order	order of the moment
   * @return	true
   */
  bool		hasMoment( int order );
  
  /**
   * @copydoc Distribution::rescale(double factor)
   */
  exponentialDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  exponentialDistribution *copy();

  /**
   * @brief Sampling from the law. Uses the mother class method.
   *
   * @return a sample from the exponential distribution
   */
  double	sample();

 public:
  /**
   * @brief Printing a representation of the law into a string
   * @author	Alain Jean-Marie
   * @return a string with a symbol for the law and its parameter
   */
  std::string toString();
  /**
   * Printing a representation of the law
   * @author	Alain Jean-Marie
   * @param out the file descriptor of the output stream
   * @param mode	representation for the output
   */
  void write( FILE *out, int mode );

};

#endif // EXPONENTIALDISTRIBUTION_H
