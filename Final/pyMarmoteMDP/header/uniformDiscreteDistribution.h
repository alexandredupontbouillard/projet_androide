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

#ifndef uniformDiscreteDistribution_H
#define uniformDiscreteDistribution_H

#include "discreteDistribution.h"

/**
 * @brief The uniform discrete distribution
 *
 */
class uniformDiscreteDistribution : public virtual discreteDistribution {

 public:
/**
 * @brief The unique constructor for the uniform discrete distribution on some interval [a..b]
 *
 * @param valInf the lower bound (a)
 * @param valSup the upper bound (b)
 */
  uniformDiscreteDistribution( int valInf, int valSup );

private:
  int _valInf; /**< @brief lower extremity of the interval */
  int _valSup; /**< @brief higher extremity of the interval */
  int _span; /**< @brief length of the interval */
  bool _isConstant; /**< @brief indicator of the fact that the interval is reduced to a point */

public:
  // accessors to specific variables
  /**
   * @brief Read accessor to the lower end of the interval
   *
   * @return the lower end of the interval
   */
  int	valInf();
  /**
   * @brief Read accessor to the upper end of the interval
   *
   * @return the upper end of the interval
   */
  int	valSup();

 public: // probabilistic member functions
  /**
    * @copydoc discreteDistribution::getProba(double value)
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
   * @copydoc Distribution::variance()
   */
  double	variance();
  /**
   * @copydoc Distribution::laplace(double)
   */
  double	laplace( double s );
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
  uniformDiscreteDistribution *rescale( double factor );
  /**
   * @copydoc Distribution::copy()
   */
  uniformDiscreteDistribution *copy();
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

#endif // uniformDiscreteDistribution_H
