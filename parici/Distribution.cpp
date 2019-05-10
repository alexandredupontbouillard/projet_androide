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

#include "header/Distribution.h"
#include "header/discreteDistribution.h"
#include "header/geometricDistribution.h"
#include <stdlib.h>
#include <iostream>
long int random() __THROW;

// constants of the class Distribution and its siblings
// MAX_RANDOM_VALUE is defined as RAND_MAX+1 so that the distribution
// normalized by it is more closely uniform on [0,1] and avoids the
// value 1. Otherwise, floor( n*u_0_1() ) does not return a uniform
// integer in the range 0..(n-1).
#define MAX_RANDOM_VALUE          0x80000000
const double Distribution::VALUE_TOLERANCE = 1.0e-8;

// calculation of the variance
double Distribution::variance()
{
  return moment(2) - ( _mean * _mean );
}

// standard implementation of iid sampling
void Distribution::iidSample( int nbSamples, double* sequence )
{
  for ( int i = 0; i < nbSamples; i++ ) {
    sequence[i] = sample();
  }
}

// the basic uniform [0,1] generator
double Distribution::u_0_1( void )
{
  return (double)random()/MAX_RANDOM_VALUE;
}

// exponential generator
double Distribution::exponential( double mean )
{
  double u = u_0_1();
  while ( u == 0.0 ) u = u_0_1();
  double res = -mean * log( u );
  /*
    if ( isinf( res ) ) fprintf( stderr, "Wow! this is oo! m=%e, u=%e, r=%f\n", 
			       mean, u, res );
  */
  return res;
}

// computing the L1 distance
double Distribution::distanceL1( Distribution *d )
{
  double res = 0.0;

  std::cout << "Ok I'm computing a distance between a " << _name
	    << " and a " << d->name() << std::endl;

  // if ( ( _name == "discreteDistribution" ) 
  //      && ( d->name() == "discreteDistribution" ) ) {
  //   discreteDistribution* d1 = dynamic_cast<discreteDistribution*>(this);
  //   discreteDistribution* d2 = dynamic_cast<discreteDistribution*>(d);

  //   res = d1->distanceL1(d2);
  // }
  // else 
  if ( hasProperty( "discreteDistribution" )
	    && d->hasProperty( "discreteDistribution" ) ) {
    // std::cout << "... mmh: both are discrete." << std::endl;
    discreteDistribution* d1 = dynamic_cast<discreteDistribution*>(this);
    discreteDistribution* d2 = dynamic_cast<discreteDistribution*>(d);

    res = d1->distanceL1(d2);
  }
  else if ( hasProperty( "discreteDistribution" )
	    && ( d->name() == "geometricDistribution" ) ) {
    discreteDistribution* d1 = dynamic_cast<discreteDistribution*>(this);
    geometricDistribution* d2 = dynamic_cast<geometricDistribution*>(d);
    for ( int i = 0; i < d1->nbVals(); i++ ) {
      res += fabs( d1->getProba(i) - d2->getProba(i) );
    }
    res += pow( d2->getRatio(), d1->nbVals() );
  }
  else if ( d->hasProperty( "discreteDistribution" )
	    && ( _name == "geometricDistribution" ) ) {
    discreteDistribution* d1 = dynamic_cast<discreteDistribution*>(d);
    geometricDistribution* d2 = dynamic_cast<geometricDistribution*>(this);
    for ( int i = 0; i < d1->nbVals(); i++ ) {
      res += fabs( d1->getProba(i) - d2->getProba(i) );
    }
    res += pow( d2->getRatio(), d1->nbVals() );
  }
  else {
    std::cout << "... not implemented." << std::endl;
  }

  return res;
}

// test of various properties.
bool Distribution::hasProperty( std::string pro )
{
  bool res = false;

  if ( pro == "discreteDistribution" ) {
    res = ( _name == "discreteDistribution" )
      || ( _name == "uniformDiscreteDistribution" )
      || ( _name == "bernoulliDistribution" )
      || ( _name == "diracDistribution" );
  }
  else if ( pro == "discrete" ) {
    res = ( _name == "discreteDistribution" )
      || ( _name == "uniformDiscreteDistribution" )
      || ( _name == "bernoulliDistribution" )
      || ( _name == "diracDistribution" )
      || ( _name == "geometricDistribution" );
  }
  else if ( pro == "continuous" ) {
    res = ( _name == "uniformDistribution" )
      || ( _name == "exponentialDistribution" );
  }
  else {
    std::cerr << "Error in Distribution::hasProperty(): property "
	      << pro << " not known. False returned in doubt." << std::endl;
  }

  return res;
}
