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

#include "multiDimHomTransition.h"
#include "../Distribution/diracDistribution.h"
#include <stdlib.h>

multiDimHomTransition::multiDimHomTransition(int nbDims, int *dimSize, double *p, double *q)
{
  _type = DISCRETE;
  _uniformizationRate = 0.0;

  _nbDims = nbDims;
  _dimSize = (int*)malloc( _nbDims*sizeof(int) );
  _p = (double*)malloc( _nbDims*sizeof(double) );
  _q = (double*)malloc( _nbDims*sizeof(double) );

  // computing the state space size and the standard self-transition probability.
  _origSize = 1;
  _r = 1.0;
  for ( int i = 0; i < _nbDims; i++ ) {
    _dimSize[i] = dimSize[i];
    _origSize *= _dimSize[i];
    _p[i] = p[i];
    _q[i] = q[i];
    _r = _r - ( p[i] + q[i] );
  }
  _destSize = _origSize;

  if ( _r < 0.0 ) {
    cerr << "Warning in multiDimHomTransition: incorrect transition probabilities. Continuing anyway."
         << endl;
  }

  // values related to offsets in the state space
  _mu = (int*)malloc( _nbDims*sizeof(int) );
  _mu[_nbDims-1] = 1;
  for ( int d = _nbDims-2; d >= 0; d-- ) {
    _mu[d] = _mu[d+1] * _dimSize[d+1];
  }

  _fromState = (int*)malloc( _nbDims*sizeof(int) );
  _toState = (int*)malloc( _nbDims*sizeof(int) );

  // buffering of transitions
  _nbMaxTrans = 1 + 2*_nbDims;
  _colNum = (int*)malloc( _nbMaxTrans*sizeof(int) );
  _colVal = (double*)malloc( _nbMaxTrans*sizeof(double) );
  _storedState = -1; // no state initially stored

}

/**
  * @brief Destructor
  *
 */
multiDimHomTransition::~multiDimHomTransition()
{
  free( _dimSize );
  free( _p );
  free( _q );
  free( _mu );
  free( _fromState );
  free( _toState );

}

/**
 * @brief Method to set the value associated with some transition. Not applicable here.
 *
 * @param i the origin state
 * @param j the destination state
 * @param val the value attached to the transition
 * @return true if the operation was successful, false otherwise (out of bounds; wrong numeric value)
 */
bool multiDimHomTransition::setEntry(int i, int j, double val)
{
  cerr << "Warning: Attempting to set entry to multiDimHomTransition transition structure. false returned."
       << endl;
  return false;
}

/**
 * @brief Method to get the value associated with some transition.
 *
 * @param i the origin state
 * @param j the destination state
 * @return the value attached to the transition (i,j)
 */
double multiDimHomTransition::getEntry(int i, int j)
{
  bool cont = true;
  double res = 0.0;

  int ires = i;
  int jres = j;
  int l1dist = 0;

  // decoding the two states dimension-by-dimension, stopping when it is clear that the
  // transition is forbidden.
  for ( int d = _nbDims-1; cont && ( d >= 0 ); d-- ) {
    _fromState[d] = ires % _dimSize[d];
    _toState[d] = jres % _dimSize[d];

    l1dist += abs( _fromState[d] - _toState[d] );
    cont = ( l1dist <= 1 );

    if ( cont ) {
      ires = ires / _dimSize[d];
      jres = jres / _dimSize[d];
    }
  }

  if ( cont ) {
    // L1-distance is 0 or 1.
    // If 0, the result depends on the location of the state wrt the boundary.
    // If 1, must find the dimension in which the move takes place
    switch( l1dist ) {
    case 0:
      res = _r;
      for ( int d = 0; d < _nbDims; d++ ) {
        if ( _fromState[d] == 0 )
          res += _q[d];
        else if ( _fromState[d] == ( _dimSize[d] - 1 ) )
          res += _p[d];
      }
      break;

    case 1:
      bool goOn = true;
      for ( int d = 0; ( d < _nbDims ) && goOn; d++ ) {
        if ( _fromState[d] != _toState[d] ) {
          if ( _toState[d] == _fromState[d] + 1 )
            res = _p[d];
          else if ( _toState[d] == _fromState[d] - 1 )
            res = _q[d];
          goOn = false;
        }
       }
      break;
    }
  }

  return res;
}

/**
 * @brief Method to get the number of non-zero entries in a transition from some state.
 * This can be seen as the number of actually possible transitions from that state.
 *
 * @param i the origin state
 * @return int the number of possible transitions
 */
int multiDimHomTransition::getNbElts(int i)
{
  int res = 1; // self jump is always possible
  int ires = i;

  // decoding state and constructing result at the same time
  for ( int d = _nbDims-1; d >= 0; d-- ) {
    _fromState[d] = ires % _dimSize[d];
    if ( _fromState[d] > 0 )
      res += 1;
    if ( _fromState[d] < ( _dimSize[d] - 1 ) )
      res += 1;

    ires = ires/_dimSize[d];
  }

  return res;
}

/**
 * @brief Method to get the @b number of the state corresponding to transition number k
 * in the list of possible transitions from some state i.
 *
 * @param i the origin state
 * @param k the index of transition from state i
 * @return the destination state corresponding to the k-th possible transition from state i
 */
int multiDimHomTransition::getCol(int i, int k)
{
  if ( i != _storedState ) {
    decodeTransitions(i);
  }

  return _colNum[k];
}

/**
 * @brief Method to get the @b value attached to transition number k
 * in the list of possible transitions from some state i.
 *
 * @param i the origin state
 * @param k the index of the transition from state i
 * @return the value attached to the k-th possible transition from state i
 */
double multiDimHomTransition::getEntryByCol(int i, int k)
{
  if ( i != _storedState ) {
    decodeTransitions(i);
  }

  return _colVal[k];
}

/**
 * @copydoc transitionStructure::getTransDistrib(int)
 */
discreteDistribution* multiDimHomTransition::getTransDistrib(int i)
{
  cerr << "Warning: method multiDimHomTransition::getTransDistrib not implemented. Dirac(0) returned." << endl;

  return new diracDistribution(0);
}

/**
 * @brief Sum of entries on some row i. Always 1.0 since this is a discrete-time
 * transition structure.
 *
 * @param i the row to be summed
 * @return 1.0
 */
double multiDimHomTransition::rowSum(int i)
{
  return 1.0;
}

/**
 * @copydoc transitionStructure::copy()
 */
multiDimHomTransition* multiDimHomTransition::copy()
{
  multiDimHomTransition* res = new multiDimHomTransition( _nbDims, _dimSize, _p, _q );

  return res;
}

/**
 * @brief Uniformizing a transition structure. Since the origin structure is already of
 * discrete-time type, a copy is returned.
 *
 * @return a discrete-time transition structure
 */
multiDimHomTransition* multiDimHomTransition::uniformize()
{
  return copy();
}

/**
 * @brief Embedding in the transition structure. Since the origin structure is already of
 * discrete-time type, a copy is returned.
 *
 * @return a discrete-time transition structure
 */
multiDimHomTransition* multiDimHomTransition::embed()
{
  return copy();
}
/**
 * @brief Computing the action of the transition structure on some measure, the measure
 * being represented as a vector of real numbers.
 * This corresponds to the multiplication vector/matrix, the row vector being interpreted
 * as a signed measure.
 * The result is placed in an array that must have been previously allocated.
 *
 * @param pi the measure to evaluate
 * @param res the resulting measure
 */
void multiDimHomTransition::evaluateMeasure(double* pi, double* res)
{
  // initialize state buffer
  decodeState( 0, _fromState );

  // initialize result
  for ( int i = 0; i < _origSize; i++ ) res[i] = 0.0;

#define TRANS( ii, jj, pro ) { res[jj] += pi[ii] * pro; }

  // enumerate all states of the parallelipiped
  for ( int i = 0; i < _origSize; i++ ) {
    // handle state number i
    double pSelf = _r;
    for ( int d = 0; d < _nbDims; d++ ) {
      if ( _fromState[d] == ( _dimSize[d] - 1 ) ) {
        pSelf += _p[d];
      }
      else {
        TRANS( i, i+_mu[d], _p[d] )
      }
      if ( _fromState[d] == 0 ) {
        pSelf += _q[d];
      }
      else {
        TRANS( i, i - _mu[d], _q[d] )
      }
    }
    TRANS( i, i, pSelf )

    // increment state
    nextState( _fromState );

  }
}

// Reimplementation of the transitionStructure generic method. I don't understand why
// I have to do that.
discreteDistribution* multiDimHomTransition::evaluateMeasure(discreteDistribution* d) {

  double *measure;
  discreteDistribution *res;

  measure = (double*)malloc( _origSize * sizeof(double) );
  evaluateMeasure( d->probas(), measure );
  res = new discreteDistribution( _origSize, d->values(), measure );

  return res;
}

/**
 * @brief Computing the action of the transition structure on some vector of values.
 * This corresponds to the multiplication matrix/vector, the column vector being interpreted
 * as a vector of values attached to the states.
 *
 * @param v the vector of values to evaluate
 * @param res the resulting vector
 */
void multiDimHomTransition::evaluateValue(double* v, double* res)
{
  // initialize state buffer
  decodeState( 0, _fromState );

  // initialize result
  for ( int i = 0; i < _destSize; i++ ) res[i] = 0.0;

#undef TRANS
#define TRANS( ii, jj, pro ) { res[ii] += v[jj] * pro; }

  // enumerate all states of the parallelipiped
  for ( int i = 0; i < _destSize; i++ ) {
    // handle state number i
    double pSelf = _r;
    for ( int d = 0; d < _nbDims; d++ ) {
      if ( _fromState[d] == ( _dimSize[d] - 1 ) ) {
        pSelf += _p[d];
      }
      else {
        TRANS( i, i+_mu[d], _p[d] )
      }
      if ( _fromState[d] == 0 ) {
        pSelf += _q[d];
      }
      else {
        TRANS( i, i - _mu[d], _q[d] )
      }
    }
    TRANS( i, i, pSelf )

    // increment state
    nextState( _fromState );

  }
}

discreteDistribution* multiDimHomTransition::getJumpDistribution()
{
  discreteDistribution* res;
  double *values = (double*)malloc( (1+2*_nbDims)*sizeof(double) );
  double *probas = (double*)malloc( (1+2*_nbDims)*sizeof(double) );

  values[0] = 0;
  probas[0] = _r;

  for ( int d = 0; d < _nbDims; d++ ) {
    values[1+2*d] = d+1;
    values[2+2*d] = -(d+1);
    probas[1+2*d] = _p[d];
    probas[2+2*d] = _q[d];
  }

  res = new discreteDistribution( 1+2*_nbDims, values, probas );

  free( values );
  free( probas );

  return res;
}

void multiDimHomTransition::write(FILE* out, string format)
{
  if ( _origSize == INFINITE_STATE_SPACE_SIZE ) {
    cerr << "Warning in multiDimHomTransition::write(): cannot write infinite chains. Ignored."
         << endl;
    return;
  }

  if ( NULL == out ) {
    return;
  }

  // Assert: chain is finite
  if ( format == "XBORNE" ) {

#define FORMAT( ii, jj, pro ) { fprintf( out, " %12e %10d", pro, jj ); }

    // initialize state buffer
    decodeState( 0, _fromState );

    // enumerate all states of the parallelipiped
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state. Note that the Rii format requires destinations to be sorted

      double pSelf = _r;
      int degree = 0;
      // scanning possible transitions to get exact degree and self transition proba
      for ( int d = 0; d <_nbDims; d++ ) {
        if ( _fromState[d] == 0 ) {
          pSelf += _q[d];
        }
        else {
          degree++;
        }
        if ( _fromState[d] == ( _dimSize[d] - 1 ) ) {
          pSelf += _p[d];
        }
        else {
          degree++;
        }
      }
      if ( pSelf > 0.0 ) {
        degree++;
      }

      fprintf( out, "%10d %9d", i, degree );

      // scanning forwards to get downwards destinations in increasing order
      for ( int d = 0; d <_nbDims; d++ ) {
        if ( _fromState[d] != 0 ) {
          FORMAT( i, i - _mu[d], _q[d] );
        }
      }
      // self loop
      if ( pSelf > 0.0 ) {
        FORMAT( i, i, pSelf );
      }
      // now scan upwards destinations in decreasing dimension number
      for ( int d = _nbDims-1; d >= 0; d-- ) {
        if ( _fromState[d] < ( _dimSize[d] - 1 ) ) {
          FORMAT( i, i + _mu[d], _p[d] );
        }
      }

      fprintf( out, "\n" );

      // increment state
      nextState( _fromState );
    }
  }
  else if ( format == "MARCA" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii, jj, pro ); }

    // initialize state buffer
    decodeState( 0, _fromState );

    // enumerate all states of the parallelipiped
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      double pSelf = _r;
      for ( int d = 0; d < _nbDims; d++ ) {
        if ( _fromState[d] == ( _dimSize[d] - 1 ) ) {
          pSelf += _p[d];
        }
        else {
          FORMAT( i, i + _mu[d], _p[d] );
        }
        if ( _fromState[d] == 0 ) {
          pSelf += _q[d];
        }
        else {
          FORMAT( i, i - _mu[d], _q[d] );
        }
      }
      FORMAT( i, i, pSelf );

      // increment state
      nextState( _fromState );
    }

  }
  else if ( format == "Ers" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii, jj, pro ); }

    fprintf( out, "discrete sparse\n" );
    fprintf( out, "%ld\n", _origSize );

    // initialize state buffer
    decodeState( 0, _fromState );

    // enumerate all states of the parallelipiped
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      double pSelf = _r;
      for ( int d = 0; d < _nbDims; d++ ) {
        if ( _fromState[d] == ( _dimSize[d] - 1 ) ) {
          pSelf += _p[d];
        }
        else {
          FORMAT( i, i + _mu[d], _p[d] );
        }
        if ( _fromState[d] == 0 ) {
          pSelf += _q[d];
        }
        else {
          FORMAT( i, i - _mu[d], _q[d] );
        }
      }
      FORMAT( i, i, pSelf );

      // increment state
      nextState( _fromState );
    }
    fprintf( out, "stop\n" );

  }
  else if ( format == "Maple" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "(%d,%d)=%e,\n", ii+1, jj+1, pro ); }
#define FORMAT2( ii, jj, pro ) { fprintf( out, "(%d,%d)=%e\n", ii+1, jj+1, pro ); }

    // initialize state buffer
    decodeState( 0, _fromState );

    // enumerate all states of the parallelipiped
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      double pSelf = _r;
      for ( int d = 0; d < _nbDims; d++ ) {
        if ( _fromState[d] == ( _dimSize[d] - 1 ) ) {
          pSelf += _p[d];
        }
        else {
          FORMAT( i, i + _mu[d], _p[d] );
        }
        if ( _fromState[d] == 0 ) {
          pSelf += _q[d];
        }
        else {
          FORMAT( i, i - _mu[d], _q[d] );
        }
      }
      // avoid comma at the end of last record
      if ( i < _origSize-1 ) {
        FORMAT( i, i, pSelf );
      }
      else {
        FORMAT2( i, i, pSelf );
      }

      // increment state
      nextState( _fromState );
    }
  }
}

/**
 * @brief utility to convert a state index into a state array
 * @author Alain Jean-Marie
 * @param index the state index
 * @param buf the state buffer to be filled
 */
void multiDimHomTransition::decodeState( int index, int* buf )
{
  int idx = index;
  for ( int d = _nbDims-1; d >= 0; d-- ) {
    buf[d] = idx % _dimSize[d];
    idx = idx / _dimSize[d];
  }
}

/**
 * @brief utility to construct the "next" state array in the lexicographical order.
 * The argument is modified.
 * @author Alain Jean-Marie
 * @param buf the state buffer to be modified
 */
void multiDimHomTransition::nextState( int* buf )
{
  buf[_nbDims-1]++;
  for ( int d =_nbDims-1; ( d >= 0 ) && ( buf[d] == _dimSize[d] ); d-- ) {
    if ( d > 0 )
  buf[d-1]++;
    buf[d] = 0;
  }
}

/**
 * @brief utility to construct transitions from a given state:
 * what states and what probabilities. The values are stored in a "cache"
 * so that multiple calls to transitionStructure::getNbElts(),
 * transitionStructure::getCol() or transitionStructure::getEntryByCol()
 * are more efficient.
 * @author Alain Jean-Marie
 * @param index the state index
 * @param buf the state buffer to be filled
 */
void multiDimHomTransition::decodeTransitions( int index )
{
  _storedState = index; // caching the current state

  _colNum[0] = index; // self jump is always possible
  _colVal[0] = _r;

  int ires = index;
  int k = 1;

  // decoding state and constructing result at the same time
  for ( int d = _nbDims-1; d >= 0; d-- ) {
    _fromState[d] = ires % _dimSize[d];
    if ( _fromState[d] > 0 ) {
      _colNum[k] = index - _mu[d];
      _colVal[k] = _q[d];
      k++;
    }
    else {
      _colVal[0] += _q[d];
    }
    if ( _fromState[d] < ( _dimSize[d] - 1 ) ) {
      _colNum[k] = index + _mu[d];
      _colVal[k] = _p[d];
      k++;
    }
    else {
      _colVal[0] += _p[d];
    }

    ires = ires/_dimSize[d];
  }

}
