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

#include "header/eventMixture.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>

using namespace std;

eventMixture::eventMixture(int size, int nbEvents, double* probas, std::string* names, int** transitions)
{
  _origSize = size;
  _destSize = _origSize;
  _type = DISCRETE;
  _uniformizationRate = 0.0;

  _nbEvents = nbEvents;
  _eventProba = (double*)malloc( _nbEvents*sizeof(double) );
  for ( int i = 0; i < _nbEvents; i++ ) {
    _eventProba[i] = probas[i];
  }
  _eventName = (std::string*)malloc( _nbEvents*sizeof(std::string));
  for ( int i = 0; i < _nbEvents; i++ ) {
    _eventName[i] = names[i];
  }
  _eventTransition = transitions;

}

// Construction from a matrix
eventMixture::eventMixture(sparseMatrix* spMat)
{
  _origSize = spMat->size();
  _destSize = _origSize;
  _type = DISCRETE;
  _uniformizationRate = 0.0;

  int* idx = (int*)malloc( _origSize*sizeof(int) );
  double* remainder = (double*)malloc( _origSize*sizeof(double) );

  _nbEvents = 0;
  double epsilon = 1e-7;
  _eventProba = (double*)malloc( _nbEvents*sizeof(double) );
  _eventTransition = (int**)malloc( _nbEvents*sizeof(int*) );

  // preparation: arrays of indexes and remainders are positioned
  // to the first element in each row.
  for ( int i = 0; i < _origSize; i++ ) {
    idx[i] = 0;
    remainder[i] = spMat->getEntryByCol(i,idx[i]);
  }

  // execute greedy event decomposition
  bool finished = false;
  do {
    // find minimal proba
    double pmin = 1.0;
    for ( int i = 0; i < _origSize; i++ ) {
      if ( remainder[i] < pmin ) pmin = remainder[i];
    }
    // record the new event
    _nbEvents++;
    _eventProba = (double*)realloc( _eventProba, _nbEvents*sizeof(double) );
    _eventProba[_nbEvents-1] = pmin;
    _eventTransition = (int**)realloc( _eventTransition, _nbEvents*sizeof(int**) );
    _eventTransition[_nbEvents-1] = (int*)malloc( _origSize*sizeof(int) );
    for ( int i = 0; i < _origSize; i++ ) {
      _eventTransition[_nbEvents-1][i] = spMat->getCol(i,idx[i]);
    }
    // remove it from remainders and advance index
    for ( int i = 0; i < _origSize; i++ ) {
      remainder[i] -= pmin;
      if ( fabs(remainder[i]) < epsilon ) {
        idx[i]++;
        if (idx[i] > spMat->getNbElts(i) ) {
          // probabilities are exhausted. Normally the algorithm terminates
          // after this loop
          finished = true;
        }
        else {
          remainder[i] = spMat->getEntryByCol(i,idx[i]);
        }
      }
    }

  } while( !finished );

  _eventName = new string[_nbEvents];
  for ( int i = 0; i < _nbEvents; i++ ) {
    std::ostringstream stmp;
    stmp << "Event" << i;
    _eventName[i] = stmp.str();
  }

  free( idx );
  free( remainder );

}

eventMixture::~eventMixture()
{
  free( _eventProba );
  free( _eventName );
}

/**
 * @brief Method to set the value associated with some transition. Not applicable here.
 *
 * @param i the origin state
 * @param j the destination state
 * @param val the value attached to the transition
 * @return false since it is not possible to set values in this structure
 */
bool eventMixture::setEntry(int i, int j, double val)
{
  cerr << "Warning: Attempting to set entry to eventMixture transition structure. false returned."
       << endl;
  return false;
}

/**
 * @copydoc transitionStructure::getEntry()
 */
double eventMixture::getEntry(int i, int j)
{
  double res = 0.0;

  cerr << "Warning: attempting to get entry to eventMixture transition structure. 0.0."
       << endl;

  return res;
}

/**
 * @copydoc transitionStructure::getNbElts()
 */
int eventMixture::getNbElts(int i)
{
  return _nbEvents;
}

/**
 * @copydoc transitionStructure::getCol()
 */
int eventMixture::getCol(int i, int k)
{
  return _eventTransition[i][k];
}

/**
 * @copydoc transitionStructure::getEntryByCol()
 */
double eventMixture::getEntryByCol(int i, int k)
{
  return _eventProba[k];
}

/**
 * @copydoc transitionStructure::getTransDistrib(int)
 */
discreteDistribution* eventMixture::getTransDistrib(int i)
{
  double* values = (double*)malloc( _nbEvents*sizeof(double) );
  for ( int i = 0; i < _nbEvents; i++ ) values[i] = i;

  discreteDistribution* res = new discreteDistribution( _nbEvents, values, _eventProba );

  return res;
}

/**
 * @brief Sum of entries on some row i. Always 1.0 since this is a discrete-time
 * transition structure.
 *
 * @param i the row to be summed
 * @return 1.0
 */
double eventMixture::rowSum(int i)
{
  return 1.0;
}

/**
 * @copydoc transitionStructure::copy()
 */
eventMixture* eventMixture::copy()
{
  // should copy the transitions!! There will be a bug when the structures are freed.
  eventMixture* res = new eventMixture( _origSize, _nbEvents, _eventProba, _eventName, _eventTransition );

  return res;
}

/**
 * @brief Uniformizing a transition structure. Since the origin structure is already of
 * discrete-time type, a copy is returned.
 *
 * @return a discrete-time transition structure
 */
eventMixture* eventMixture::uniformize()
{
  return copy();
}

/**
 * @brief Embedding in a transition structure. Since the origin structure is already of
 * discrete-time type, a copy is returned.
 *
 * @return a discrete-time transition structure
 */
eventMixture* eventMixture::embed()
{
  return copy();
}

void eventMixture::evaluateMeasure(double *d, double *res)
{
  // initialize result
  for ( int j = 0; j < _destSize; j++ ) res[j] = 0.0;

  // scan all events and transitions in each event
  for ( int k = 0; k < _nbEvents; k++ ) {
    for ( int i = 0; i < _origSize; i++ ) {
      res[ _eventTransition[k][i] ] += d[i] * _eventProba[k];
    }
  }
}

discreteDistribution* eventMixture::evaluateMeasure(discreteDistribution* d) {

  double *measure;
  discreteDistribution *res;

  measure = (double*)malloc( _destSize * sizeof(double) );
  evaluateMeasure( d->probas(), measure );
  res = new discreteDistribution( _destSize, d->values(), measure );

  return res;
}

/**
 * @copydoc transitionStructure::evaluateValue(double*,double*)
 */
void eventMixture::evaluateValue(double *v, double *res)
{
  // initialize result
  for ( int i = 0; i < _origSize; i++ ) res[i] = 0.0;

  // scan all events and transitions in each event
  for ( int k = 0; k < _nbEvents; k++ ) {
    for ( int i = 0; i < _origSize; i++ ) {
      res[i] += v[ _eventTransition[k][i] ] * _eventProba[k];
    }
  }
}

double eventMixture::evaluateValueState(double *v, int stateIndex)
{
  double res = 0.0;

  // scan all events and transitions in each event
  for ( int k = 0; k < _nbEvents; k++ ) {
    res += v[ _eventTransition[k][stateIndex] ] * _eventProba[k];
  }

  return res;
}

/**
 * @brief Output method for the transition structure. Supported formats are:
 * XBORNE, MARCA, Ers, Maple.
 *
 * @param out the file descriptor to which the structure should be written.
 * @param format the format/language to be used
 */
void eventMixture::write(FILE* out, string format)
{

  if ( NULL == out ) {
    return;
  }

  // arrays necessary for consolidating transition probabilities
  int* destinations = (int*)malloc( _nbEvents*sizeof(int) );
  double* values = (double*)malloc( _nbEvents*sizeof(double) );

  if ( format == "XBORNE" ) {

#define FORMAT( ii, jj, pro ) { fprintf( out, " %12e %10d", pro, jj ); }

    // enumerate all states
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      int nbTrans = consolidate( i, destinations, values );
      fprintf( out, "%10d %10d", i, nbTrans );
      for ( int d = 0; d < nbTrans; d++ ) {
        FORMAT( i, destinations[d], values[d] );
      }
      fprintf( out, "\n" );
    }
  }

  else if ( format == "MARCA" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii, jj, pro ); }

    // enumerate all states
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      int nbTrans = consolidate( i, destinations, values );
      for ( int d = 0; d < nbTrans; d++ ) {
        FORMAT( i, destinations[d], values[d] );
      }
    }
  }
  else if ( format == "Ers" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii, jj, pro ); }

    // fprintf( out, "discrete sparse\n" );
    // fprintf( out, "%ld\n", _origSize );

    // enumerate all states
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      // enumerate events
      for ( int d = 0; d < _nbEvents; d++ ) {
        FORMAT( i, _eventTransition[d][i], _eventProba[d] );
      }
    }
    // fprintf( out, "stop\n" );
  }
  else if ( format == "Maple" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "(%d,%d)=%e,\n", ii+1, jj+1, pro ); }
#define FORMAT2( ii, jj, pro ) { fprintf( out, "(%d,%d)=%e\n", ii+1, jj+1, pro ); }

    // enumerate all states
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      int nbTrans = consolidate( i, destinations, values );
      if ( i < _origSize-1 ) {
        for ( int d = 0; d < nbTrans; d++ ) {
          FORMAT( i, destinations[d], values[d] );
        }
      }
      else {
        for ( int d = 0; d < nbTrans-1; d++ ) {
          FORMAT( i, destinations[d], values[d] );
        }
        // avoid comma at the end of last record
        FORMAT2( i, destinations[nbTrans-1], values[nbTrans-1] );
      }
    }
  }
  else if ( format == "PSI3" )
  {
    string fileName = "param.yaml";
    FILE *out;

    if ( NULL == ( out = fopen( fileName.c_str(), "w" ) ) ) {
      cerr << "Error in eventMixture::write(): "
           << "could not open output file '" << fileName << "'. Ignored" << endl;
    }
    else {
      fprintf( out, "%-24s%s\n", "Version:", "1.0" );
      fprintf( out, "\n" );
      fprintf( out, "%-24s%s\n", "MyLib:", "[./mylib/stopfct]" );
      fprintf( out, "\n" );
      fprintf( out, "%-24s%s\n", "Method:", "simpleforward" );
      fprintf( out, "\n" );
      fprintf( out, "%-24s%d\n", "SampleNumber:", 1000 );
      fprintf( out, "\n" );
      fprintf( out, "%s\n", "# Stop criteria of the simulation algorithm" );
      fprintf( out, "%-24s%s\n", "StopFct:", "Default$stop" );
      fprintf( out, "\n" );
      fprintf( out, "%s\n", "#Random generator seed. \"~\" for random seed" );
      fprintf( out, "%-24s%s\n", "Seed:", "~" );
      fprintf( out, "\n" );
      fprintf( out, "%s\n", "#Initial states of diffents queues of the simulation algorithm" );
      fprintf( out, "%s\n", "#Default value is chosen randomly (\" ~ \")" );
      fprintf( out, "%-24s%s\n", "InitialState :", "[0]" );
    }
  }
  else {
    cerr << "Warning: in eventMixture::write(): format '" << format
         << "' not supported. Ignored." << endl;
  }

  free( destinations );
  free( values );
}

// method to consolidate transitions from some state
int eventMixture::consolidate(int i, int* destinations, double* values)
{
  int nbTrans = 0;

  for ( int e = 0; e < _nbEvents; e++ ) values[e] = 0.0;

  // straighforward insertion. No sorting, no shifting.
  for ( int e = 0; e < _nbEvents; e++ ) {
    int dest = _eventTransition[i][e];
    bool isNew = true;
    int location = nbTrans;
    for ( int i = 0; isNew && ( i < nbTrans ); i++ ) {
      isNew = ( dest == destinations[i] );
      location = i;
    }
    if ( isNew ) {
      destinations[location] = dest;
      values[location] = _eventProba[e];
    }
    else {
      values[location] += _eventProba[e];
    }
  }

  return nbTrans;
}
