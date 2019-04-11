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

#ifndef geometricDistribution_H
#define geometricDistribution_H

#include "Distribution.h"

/**
 * @brief The geometric distribution with starting value 0. The parameter "p" is called "ratio".
 * The Geometric distribution is discrete but does not inherit from discreteDistribution
 * because its range is infinite.
 * @author Alain Jean-Marie
 *
 */
class geometricDistribution : public virtual Distribution {

public:
  // constructors
/**
 * @brief Unique constructor for the class, from its "ratio".
 * @author Alain Jean-Marie
 * @param p the probability of being larger than 0
 */
  geometricDistribution( double p );

private:
  // specific variables
  double _p; /**< Value of the parameter, which is the probability of being larger than 0 */

public:
  // accessors
  /**
   * @brief Function to obtain the probability of a specific value k
   *
   * @param k the value at which the probability should be computed
   * @return the probability that the random variable is k
   */
  double 	getProba(double k );
  /**
   * @brief Function to obtain the parameter (or ratio) of the distribution. Redundant with p() but
   * defined to be more explicit.
   *
   * @return the value of the ratio
   */
  double 	p() { return _p; }
  /**
    * @brief Function to obtain the parameter (or ratio) of the distribution. Redundant with getRatio() but
    * defined according to the coding convention.
    *
    * @return double
    */
  double 	getRatio() { return _p; }

public:
  // probabilistic member functions
  /**
   * @brief Function to obtain the mean (expectation). Its value is 1/(1-_p) @see geometricDistribution::_p
   *
   * @return the mathematical expectation of the distribution
   */
  double	mean();
  /**
   * @brief
   *
   * @return double
   */
  double	rate(); // inverse of the mean
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
   * @brief
   *
   * @param x
   * @return double
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
   * @brief Rescaling the distribution. Geometric distributions cannot be
   * rescaled. A copy is returned. A warning is issued if the factor is not 1.0.
   * @param factor the factor by which to rescale
   * @return the rescaled distribution.
   */
  geometricDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  geometricDistribution *copy();

  /**
   * @brief Sampling from the distribution. The method uses the fact that
   * the integer part of an exponential random variable is a geometric random
   * variable.
   *
   * @return a sample from the geometric distribution
   */
  double	sample();

 public:
  /**
   * @brief
   *
   * @return std::string
   */
  std::string toString();
  /**
   * @brief
   *
   * @param out
   * @param mode
   */
  void write( FILE *out, int mode );

};

#endif // geometricDistribution_H
