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

#ifndef diracDistribution_H
#define diracDistribution_H

#include "discreteDistribution.h"

/**
 * @brief The Dirac distribution concentrated at some point

 *
 */
class diracDistribution : public virtual discreteDistribution {

public:
/**
 * @brief Unique constructor for the Dirac distribution from its value.
 *
 * @param val the value at which the distribution is concentrated
 */
  diracDistribution( double val ); // creation from a value

private:
  // specific variables
  double _value; /**< the value characterizing the distribution */

public:
  // accessors to specific variables
  /**
   * @brief Read accessor to the value of the Dirac distribution
   *
   * @return the value of the distribution
   */
  double	value(int) { return _value; };

public: 
  // probabilistic member functions

  /**
   * @copydoc discreteDistribution::mean()
   */
  double getProba(double value);

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
   * @brief Variance of the Dirac distribution. Reimplemented to return
   * always 0.
   */
  double	variance();
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
   * @copydoc Distribution::rescale(double)
   */
  diracDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  diracDistribution *copy();

  /**
   * @copydoc Distribution::sample()
   */
  double	sample();
  /**
   * @brief Sampling of i.i.d. values in a table. Reimplemented in order
   * to avoid useless function calls.
   *
   * @copydetails Distribution::iidSample(int n, double* s)
   */
  void		iidSample( int n, double* s );

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

#endif // diracDistribution_H
