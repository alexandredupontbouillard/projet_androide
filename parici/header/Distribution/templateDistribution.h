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

#ifndef templateDistribution_H
#define templateDistribution_H

#include "Distribution.h"

/**
 * @brief The general template distribution to be instantiated
 *
 */
class templateDistribution : public virtual Distribution {

public:
  /**
   * @brief Constructor for some distribution of some type from things and other.
   * The mean is calculated at creation (or not).
   * @author	Alain Jean-Marie
   * @param x un truc
   * @param y un machin
   * @param z une chose
   * @return	an object of type templateDistribution
   */
  templateDistribution( int x, double* y, double* z ); // creation from something

private:
  // private variables specific to the distribution
  int	  _nbVals; /**< typicaly: number of values */
  double* _values; /**< typically: the values */
  double* _probas; /**< typically: probabilities attached to values */

public:
  // accessors to specific variables

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
   * @copydoc Distribution::variance()
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
   * @copydoc Distribution::ccdf()
   */
  double	ccdf( double x );
  /**
   * @copydoc Distribution::hasMoment(int)
   */
  bool		hasMoment( int order );
  
  /**
   * @copydoc Distribution::rescale(double)
   */
  templateDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  templateDistribution *copy();
  /**
   * @copydoc Distribution::sample()
   */
  double	sample();
  /**
   * @copydoc Distribution::iidSample()
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

#endif // templateDistribution_H
