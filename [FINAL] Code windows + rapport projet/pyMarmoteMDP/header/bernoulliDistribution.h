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

#ifndef bernoulliDistribution_H
#define bernoulliDistribution_H

#include "discreteDistribution.h"

/**
 * @brief The Bernoulli distribution with two values
 *
 */
class bernoulliDistribution : public virtual discreteDistribution {

public:
  /**
   * @brief Unique onstructor for a Bernoulli distribution from the probability
   * that it is equal to 1.
   * @author	Alain Jean-Marie
   * @param	val the value P( X = 1 )
   * @return	an object of type bernoulliDistribution
   */
  bernoulliDistribution( double ); // creation from a value

private:
  // private variables specific to the distribution
  double _proba;  /**< proba of the value 1 */

public:
  // accessors to specific variables
  /**
   * @brief Accessor to the parameter of the distribution. Redundant with the
   * standard accessor proba() but more explicit.
   *
   * @return the parameter
   */
  double	getParameter() { return _proba; }
  /**
   * @brief Accessor to the parameter of the distribution. Redundant with the
   * accessor getParameter() but conform to the coding standard.
   *
   * @return the parameter
   */
  double	proba() { return _proba; }

public: 
  // probabilistic member functions
  /**
   * @copydoc Distribution::mean()
   */
  double	mean();
  /**
   * @copydoc Distribution::rate()
   */
  double	rate(); 
  /**
   * @copydoc Distribution::moment(int)
   */
  double	moment( int order );
  /**
   * @copydoc Distribution::laplace(double)
   */
  double	laplace( double s ); // Laplace transform at real points
  /**
   * @copydoc Distribution::dLaplace(double)
   */
  double	dLaplace( double s ); // derivative of the Laplace transform
  /**
   * @copydoc Distribution::cdf()
   */
  double	cdf( double x );
  /**
   * @copydoc Distribution::hasMoment(int)
   */
  bool		hasMoment( int order );
  
  /**
   * @brief Rescaling the distribution. Bernoulli distributions cannot be
   * rescaled. Ac copy is returned and an error message is issued if the
   * factor is not 1.0.
   *
   * @copydetails Distribution::rescale(double)
   */
  bernoulliDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  bernoulliDistribution *copy();

  /**
   * @copydoc Distribution::sample()
   */
  double	sample();

 public:
  /**
   * @copydoc Distribution::toString()
   */
  std::string toString();
  /**
   * @copydoc Distribution::write(FILE*,int)
   */
  void write( FILE *out, int mode );

};

#endif // bernoulliDistribution_H
