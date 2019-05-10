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

#ifndef poissonDistribution_H
#define poissonDistribution_H

#include "Distribution.h"
#ifdef WITH_R
#include <RInside.h>
#endif

/**
 * @brief The Poisson distribution. The parameter is called "lambda".
 * The Poisson distribution is discrete but does not inherit from discreteDistribution
 * because its range is infinite.
 * @author Alain Jean-Marie
 *
 */
class poissonDistribution : public virtual Distribution {

public:
  /**
   * @brief Constructor for a Poisson distribution from its "lambda" parameter.
   * The mean is calculated at creation.
   * @author	Alain Jean-Marie
   * @param lambda a non-negative real number
   * @return	an object of type poissonDistribution
   */
  poissonDistribution( double lambda );

  /**
   * @brief Destructor for a Poisson distribution.
   * @author	Alain Jean-Marie
   */
  ~poissonDistribution();

private:
  // private variables specific to the distribution
  double  _lambda; /**< the "rate" of the Poisson distribution, also its mean. */

#ifdef WITH_R
  static RInside* _Rmotor; /**< @brief R execution engine; there is at most one; may not be created at all */
#else
  typedef void* RInside;
  static RInside* _Rmotor;
#endif

private:
  RInside* Rmotor();


public:
  // accessors to specific variables
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
  double 	lambda() { return _lambda; }

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
  poissonDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  poissonDistribution *copy();
  /**
   * @brief Sample from the Poisson distibution. Uses the "R" package.
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

#endif // poissonDistribution_H
