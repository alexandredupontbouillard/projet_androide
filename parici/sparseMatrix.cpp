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

#include <stdlib.h>
#include "header/sparseMatrix.h"
#include "header/diracDistribution.h"

#include <sstream>
#include <string>
#include <iomanip> 


/**
 * Standard constructor for square sparse matrices. The internal structures
 * are initialized. Type is set to UNKNOWN. The result is the null matrix.
 */
sparseMatrix::sparseMatrix(int size)
{
  _origSize = (long int)size;
  _destSize = _origSize;
  _debug = false;
  _type = UNKNOWN;

  // allocation
  _nbElts = (int*) malloc( (size_t)_origSize * sizeof(int) );
  _elts = (int**) malloc( (size_t)_origSize * sizeof(int*) );
  _vals = (double**) malloc( (size_t)_origSize * sizeof(double*) );

  _sccs.first = (std::vector<SCC> *)NULL;
  _sccs.second = (sparseMatrix *)NULL;
  _reverted = (sparseMatrix *)NULL;

  // initialization of contents
  if ( ( (int*)NULL != _nbElts )
       && ( (int**)NULL != _elts )
       && ( (double**)NULL != _vals ) ) {
    for ( int i=0; i<_origSize; i++ ) {
      _nbElts[i] = 0;
      _elts[i] = (int*)NULL;
      _vals[i] = (double*)NULL;
    }
  } 
  else {
    // abort
    if ( (int*)NULL != _nbElts ) free( _nbElts );
    if ( (int**)NULL != _elts ) free( _elts );
    if ( (double**)NULL != _vals ) free( _vals );
    // should throw some exception here
  }

}

// Constructor for non-square sparse matrices.
sparseMatrix::sparseMatrix(int rowSize, int colSize)
{
  _origSize = (long int)rowSize;
  _destSize = (long int)colSize;
  _debug = false;
  _type = UNKNOWN;

  // allocation
  _nbElts = (int*) malloc( (size_t)_origSize * sizeof(int) );
  _elts = (int**) malloc( (size_t)_origSize * sizeof(int*) );
  _vals = (double**) malloc( (size_t)_origSize * sizeof(double*) );

  _sccs.first = (std::vector<SCC> *)NULL;
  _sccs.second = (sparseMatrix *)NULL;
  _reverted = (sparseMatrix *)NULL;

  // initialization of contents
  if ( ( (int*)NULL != _nbElts )
       && ( (int**)NULL != _elts )
       && ( (double**)NULL != _vals ) ) {
    for ( int i=0; i<_origSize; i++ ) {
      _nbElts[i] = 0;
      _elts[i] = (int*)NULL;
      _vals[i] = (double*)NULL;
    }
  }
  else {
    // abort
    if ( (int*)NULL != _nbElts ) free( _nbElts );
    if ( (int**)NULL != _elts ) free( _elts );
    if ( (double**)NULL != _vals ) free( _vals );
    // should throw some exception here
  }

}

/**
 * @brief Standard destructor
 * @author Alain Jean-Marie
 */
sparseMatrix::~sparseMatrix()
{
  if ( (int*)_nbElts != NULL ) {
    free( _nbElts );
  }
  if ( (int**)_elts != NULL ) {
    for ( int i = 0; i < _origSize; i++ ) free( _elts[i] );
    free( _elts );
  }
  if ( (int**)_vals != NULL ) {
    for ( int i = 0; i < _origSize; i++ ) free( _vals[i] );
    free( _vals );
  }
  delete _sccs.first;
  _sccs.first = (vector<SCC> *) NULL;
  delete _sccs.second;
  _sccs.second = (sparseMatrix *) NULL;
  delete _reverted;
  _reverted = (sparseMatrix *) NULL;
}

// Inserting an element in the matrix.
bool sparseMatrix::setEntry(int row, int col, double val)
{
  //clearComputedSCCs();
  // fprintf( stderr, "R=%d C=%d V=%f\n", row, col, val );
  // checking ranges
  if ( ( 0 <= row ) && ( row < _origSize )
       && ( 0 <= col ) && ( col < _destSize ) ) {
    _nbElts[row]++;
    _elts[row] = (int*)realloc( _elts[row],
             (size_t)_nbElts[row]*sizeof(int) );
    _vals[row] = (double*)realloc( _vals[row],
            (size_t)_nbElts[row]*sizeof(double) );

    // checking that allocation was successful
    if ( ( (int*)NULL == _elts[row] )
     || ( (double*)NULL == _vals[row] ) ) {
      return false;
    }
    else {
      _elts[row][_nbElts[row]-1] = col;
      _vals[row][_nbElts[row]-1] = val;
      return true;
    }
  }
  else {
    printf("Warning: setting value to (%d,%d) in matrix of size (%ld,%ld).",
       row, col, _origSize, _destSize );
    printf(" Ignored.\n");
    return false;
  }
}

// adding a value to the matrix
bool sparseMatrix::addToEntry(int row, int col, double val)
{
  //clearComputedSCCs();
  int i;

  // checking ranges
  if ( ( 0 <= row ) && ( row < _origSize )
       && ( 0 <= col ) && ( col < _destSize ) ) {
    // checking if the element is there already
    for ( i=0; 
          ( i<_nbElts[row] ) && ( _elts[row][i]!= col );
          i++ );
    if ( i<_nbElts[row] ) {
      // found it.
      _vals[row][i] += val;
      return true;
    }
    else { 
      /* not there. add element. I replicate the code of "setEntry" */
      /* for the sake of efficiency.					*/
      _nbElts[row]++;
      _elts[row] = (int*)realloc( _elts[row],
               (size_t)_nbElts[row]*sizeof(int) );
      _vals[row] = (double*)realloc( _vals[row],
              (size_t)_nbElts[row]*sizeof(double) );
      
      if ( ( (int*)NULL == _elts[row] )
       || ( (double*)NULL == _vals[row] ) ) {
	return false;
      }
      else {
    _elts[row][_nbElts[row]-1] = col;
    _vals[row][_nbElts[row]-1] = val;
	return true;
      }
    }
  }
  else {
    printf("Warning: setting value to (%d,%d) in matrix of size (%ld,%ld).",
       row, col, _origSize, _destSize );
    printf(" Ignored.\n");
    return false;
  }
}

// retrieving the entry of a matrix
double sparseMatrix::getEntry(int row, int col)
{
  double	res = 0.0;

  // the entire array has to be scanned because column index may be replicated
  for( int i=0; i<_nbElts[row]; i++ ) {
    if ( _elts[row][i] == col ) {
      res += _vals[row][i];
    }
  }

  return res;
}

// get the number of entries in row i
int sparseMatrix::getNbElts(int row)
{
  int res = 0;

  if ( ( row < 0 ) || ( row > _origSize ) ) {
    fprintf( stderr, "Warning in sparseMatrix::getNbElts:" );
    fprintf( stderr, " index %d out of bounds. 0 returned.\n",
	     row );
  }
  else {
    res = _nbElts[row];
  }

  return res;
}

// get the jth column for row i
int sparseMatrix::getCol(int row, int numCol)
{
  int res = -1;

  if ( ( numCol < 0 ) || ( numCol > _nbElts[row] ) ) {
    fprintf( stderr, "Warning in sparseMatrix: index %d out of bounds. -1 returned.\n",
	     numCol );
  }
  else {
    res = _elts[row][numCol];
  }

  return res;
}

// get the value in the jth column for row i
double sparseMatrix::getEntryByCol(int row, int numCol)
{
  double res = 0.0;

  if ( ( numCol < 0 ) || ( numCol > _nbElts[row] ) ) {
    fprintf( stderr, "Warning in sparseMatrix::getEntryBycol(): index %d out of bounds. 0 returned.\n",
	     numCol );
  }
  else {
    res = _vals[row][numCol];
  }

  return res;
}

// transitions from some state
discreteDistribution* sparseMatrix::getTransDistrib(int state)
{
  double* next;
  double* probas;
  discreteDistribution *res;

  switch ( _type ) {
  case DISCRETE:
    // for a discrete-time structure, it is not supposed to happen
    // that the number of elements is 0.
    next = (double*) malloc( (size_t)_nbElts[state] * sizeof(double) );
    probas = (double*) malloc( (size_t)_nbElts[state] * sizeof(double) );

    for ( int j = 0; j < _nbElts[state]; j++ ) {
      next[j] = _elts[state][j];
      probas[j] = _vals[state][j];
    }
    res = new discreteDistribution( _nbElts[state], next, probas );

    // cleanup temporary arrays: they have been copied at creation of discreteDistribution
    free( next );
    free( probas );

    break;

  case CONTINUOUS:
    // for a continuous-time structure, it may happen that the number of
    // elements is 0. Since there is no discrete distribution with 0 values,
    // we force the distribution to be Dirac at i.
    if ( _nbElts[state] == 0 ) {
      res = new diracDistribution( (double)state );
    }
    else {
      // assumption that the diagonal entry is unique in the list...
      next = (double*) malloc( (size_t)(_nbElts[state]-1) * sizeof(double) );
      probas = (double*) malloc( (size_t)(_nbElts[state]-1) * sizeof(double) );

      double totRate = -getEntry(state,state);
      int idx = 0;
      for ( int j = 0; j < _nbElts[state]; j++ ) {
        if ( _elts[state][j] != state ) {
          next[idx] = _elts[state][j];
          probas[idx] = _vals[state][j] / totRate;
          idx++;
        }
      }

      res = new discreteDistribution( _nbElts[state]-1, next, probas );

      // cleanup temporary arrays: they have been copied at creation of discreteDistribution
      free( next );
      free( probas );
    }
    break;

  case UNKNOWN:
    fprintf( stderr, "Warning in getTransDistrib(): unknown type for sparseMatrix. Dirac(0) returned.\n");
    res = new diracDistribution(0);
  }

  return res;
}

// summing a row
double sparseMatrix::rowSum(int row)
{
  double sum = 0.0;

  for( int i=0; i<_nbElts[row]; i++) {
    sum += _vals[row][i];
  }

  return sum;
}

/* Linear algebra... */
// multiplication vector/matrix.
void sparseMatrix::evaluateMeasure(double* m, double *res)
{
  int		k;

  // (re) initialization
  for ( int i=0; i<_origSize; i++ )
    res[i] = 0.0;

  for ( int i=0; i<_origSize; i++ ) {
    for ( int j=0; j<_nbElts[i]; j++ ) {
      k = _elts[i][j];
      res[k] += m[i] * _vals[i][j];
    }
  }
}

// Reimplementation of the transitionStructure generic method. I don't understand why
// I have to do that.
discreteDistribution* sparseMatrix::evaluateMeasure(discreteDistribution* d) {

  double *measure;
  discreteDistribution *res;

  measure = (double*)malloc( _origSize * sizeof(double) );
  evaluateMeasure( d->probas(), measure );
  res = new discreteDistribution( _origSize, d->values(), measure );

  return res;
}

// multiplication matrix/vector
void sparseMatrix::evaluateValue(double* v, double *res)
{
  for ( int i=0; i<_destSize; i++ ) {
    res[i] = 0.0;
    for ( int j=0; j<_nbElts[i]; j++ ) {
      res[i] += _vals[i][j] * v[ _elts[i][j] ];
    }
  }
}

// multiplication matrix row/vector for a given state.
double sparseMatrix::evaluateValueState(double* v, int stateIndex)
{
	//evaluate only for a givin state and return results.
	double res = 0.0;
	for ( int j=0; j<_nbElts[stateIndex]; j++ )
	{
	      res += _vals[stateIndex][j] * v[ _elts[stateIndex][j] ];
	}
	return res;
}

// Standard copy.
sparseMatrix* sparseMatrix::copy()
{
  sparseMatrix* res;

  res = new sparseMatrix( _origSize, _destSize );
  res->setType( _type );

  for ( int i=0; i < _origSize; i++ ) {
    for ( int j=0; j < _nbElts[i]; j++ ) {
      res->setEntry( i, j, _vals[i][j] );
    }
  }

  return res;

}

// Uniformization
sparseMatrix* sparseMatrix::uniformize()
{
  sparseMatrix* res;
  double lambda;
  double nu = 0.0;

  if ( _type == UNKNOWN ) {
    fprintf( stderr, "Error: cannot uniformize matrix of unknown type. Copy of the original returned.\n" );
    res = copy();
  }
  else if ( _type != CONTINUOUS ) {
    fprintf( stderr, "Error: cannot uniformize discrete-time matrix. Copy of the original returned.\n" );
    res = copy();
  }
  else if ( _origSize != _destSize ) {
    fprintf( stderr, "Error: cannot uniformize non-square matrix. Null returned.\n" );
    res = (sparseMatrix*)NULL;
  }
  else if ( _origSize == INFINITE_STATE_SPACE_SIZE ) {
    fprintf( stderr, "Error: cannot uniformize sparse matrices with infinite state spaces. Null returned.\n" );
    res = (sparseMatrix*)NULL;
  }
  else {
    // assert: continuous-time generator with finite state space.
    // Computing the uniformization factor
    for ( int i=0; i<_origSize; i++ ) {
      lambda = - getEntry( i, i );
      nu = ( lambda > nu ? lambda : nu );
    }
    
    sparseMatrix* uGen = new sparseMatrix( _origSize );
    uGen->setType( DISCRETE );
    uGen->setUniformizationRate( nu );

    for ( int i=0; i < _origSize; i++ ) {
      if ( _nbElts[i] == 0 ) { // this is an absorbing state
	uGen->setEntry( i, i, 1.0 );
      }

      for ( int j=0; j < _nbElts[i]; j++ ) {
        if ( _elts[i][j] != i )
          uGen->setEntry( i, _elts[i][j], _vals[i][j] / nu );
        else
          uGen->setEntry( i, i, 1.0 + _vals[i][j] / nu );
      }
    }

    if ( _debug ) {
      fprintf( stdout, "Uniformized with factor nu = %f:\n", nu );
      uGen->write( stdout, "Ers" );
    }

    res = uGen;
  }

  return res;
}

// Embedding
sparseMatrix* sparseMatrix::embed()
{
  sparseMatrix* res;

  if ( _type == UNKNOWN ) {
    fprintf( stderr, "Error: cannot embed matrix of unknown type. Copy of the original returned.\n" );
    res = copy();
  }
  else if ( _type != CONTINUOUS ) {
    fprintf( stderr, "Error: cannot embed a stochastic matrix. Copy of the original returned.\n" );
    res = copy();
  }
  else if ( _origSize != _destSize ) {
    fprintf( stderr, "Error: cannot embed in non-square matrix. Null returned.\n" );
    res = (sparseMatrix*)NULL;
  }
  else if ( _origSize == INFINITE_STATE_SPACE_SIZE ) {
    fprintf( stderr, "Error: cannot embed sparse matrices with infinite state spaces. Null returned.\n" );
    res = (sparseMatrix*)NULL;
  }
  else {
    // assert: continuous-time generator with finite state space.

    sparseMatrix* eGen = new sparseMatrix( _origSize );
    eGen->setType( DISCRETE );

    for ( int i=0; i < _origSize; i++ ) {
      double nu = - getEntry( i, i );
      if ( nu == 0.0 ) {
	fprintf( stderr, "Warning in embed(): state %d is absorbing.", i );
	fprintf( stderr, "Embedding will not work correctly.\n" );
	eGen->setEntry( i, i, 1.0 );
      }
      else {
	for ( int j=0; j < _nbElts[i]; j++ ) {
	  if ( _elts[i][j] != i )
	    eGen->setEntry( i, _elts[i][j], _vals[i][j] / nu );
	}
      }
    }

    if ( _debug ) {
      fprintf( stdout, "Embedding at jump times\n" );
      eGen->write( stdout, "Ers" );
    }

    res = eGen;
  }

  return res;
}

// Diagnostic method
void sparseMatrix::diagnose( FILE* out )
{

  fprintf( out, "Diagnostic for sparseMatrix structure:\n" );
  fprintf( out, "- generator type:        " );
  switch( _type ) {
  case DISCRETE:
    fprintf( out, "discrete"); break;
  case CONTINUOUS:
    fprintf( out, "continuous"); break;
  case UNKNOWN:
    fprintf( out, "unknown"); break;
  default:
    fprintf( out, "undefined");
  }
  fprintf( out, "\n" );
  fprintf( out, "- number of origin states:      %ld\n", _origSize );
  fprintf( out, "- number of destination states: %ld\n", _destSize );

  int* colhit = (int*) calloc( _destSize, sizeof(int) );
  
  // scanning transitions
  int nbVoid = 0;
  int nbTrans = 0;
  int min_oDeg = _destSize;
  int max_oDeg = 0;
  int min_iDeg = _origSize;
  int max_iDeg = 0;
  double minSum = 1e38;
  double maxSum = -1e38;
  double minElt = 1e38;
  double maxElt = -1e38;
  double mismatch = 0.0;
  for ( int i = 0; i < _origSize; i++ ) {
    if ( _nbElts[i] == 0 ) {
      nbVoid++;
    }
    else {
      int deg = _nbElts[i];
      nbTrans+= deg;
      if ( deg > max_oDeg ) max_oDeg = deg;
      if ( deg < min_oDeg ) min_oDeg = deg;
      double rs = rowSum(i);
      if ( rs < minSum ) minSum = rs;
      if ( rs > maxSum ) maxSum = rs;
      double tot = 0.0;
      for ( int j = 0; j < _nbElts[i]; j++ ) {
        double val = _vals[i][j];
	int col = _elts[i][j];
        tot += val;
        if ( val < minElt ) minElt = val;
        if ( val > maxElt ) maxElt = val;
	if ( col != i ) colhit[col]++;
      }
      double mis = fabs( tot - rs );
      if ( mis > mismatch) mismatch = mis;
    }
  }

  int nbVoidCol = 0;
  for ( int i = 0; i < _destSize; i++ ) {
    if ( colhit[i] == 0 ) nbVoidCol++;
    if ( colhit[i] > max_iDeg ) max_iDeg = colhit[i];
    if ( colhit[i] < min_iDeg ) min_iDeg = colhit[i];
  }
  
  // output diagnostic
  fprintf( out, "- number of transitions: %d\n", nbTrans );
  fprintf( out, "- number of empty rows:  %d\n", nbVoid );
  fprintf( out, "- maximum outdegree:     %d\n", max_oDeg );
  fprintf( out, "- minimum outdegree:     %d\n", min_oDeg );
  fprintf( out, "- maximum indegree:      %d\n", max_iDeg );
  fprintf( out, "- minimum indegree:      %d\n", min_iDeg );
  fprintf( out, "- maximum value:         %12.6e\n", maxElt );
  fprintf( out, "- minimum value:         %12.6e\n", minElt );
  fprintf( out, "- maximum row sum:       %12.6e\n", maxSum );
  fprintf( out, "- minimum row sum:       %12.6e\n", minSum );
  fprintf( out, "- row sum mismatch:      %12.6f\n", mismatch );

  // cleanup
  free( colhit );
}

// General output method
void sparseMatrix::write( FILE* out, std::string format )
{

  if ( NULL == out ) {
    return;
  }

  // Arrays necessary for consolidating transition probabilities.
  // Consolidating is necessary because many formats do not allow
  // multiple entry declaration with the convention that they add up.
  int maxNbElts = 0;
  for ( int i = 0; i < _origSize; i++ )
    if ( _nbElts[i] > maxNbElts ) maxNbElts = _nbElts[i];
  int* destinations = (int*)malloc( maxNbElts*sizeof(int) );
  double* values = (double*)malloc( maxNbElts*sizeof(double) );

  if (format == "XBORNE") {

    fprintf( out, "%s", toString("XBORNE").c_str() );

    // need to create modelName.sz (and modelName.cd)
  }
  else if ( ( format == "MARCA" ) || ( format == "MatrixMarket-sparse" ) ) {
    // both MARCA and Matrix Market are "1-based"

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii+1, jj+1, pro ); }

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

    for ( int i=0; i < _origSize; i++ ) {
      for ( int j=0; j < _nbElts[i]; j++ ) {
        FORMAT( i, _elts[i][j], _vals[i][j] );
      }
    }
  }
  else if ( format == "Maple" ) {
    // Maple matrices are "1-based".

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "(%d,%d)=%12e,\n", ii+1, jj+1, pro ); }
#define FORMAT2( ii, jj, pro ) { fprintf( out, "(%d,%d)=%12e\n", ii+1, jj+1, pro ); }

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
  else if ( format == "MatrixMarket-full" ) {
    // format is "column-major" (!)
    for ( int j=0; j<_destSize; j++ ) {
      for ( int i=0; i<_origSize; i++ ) {
        double val = getEntry( i, j );
        fprintf( out, " %12e", val );
      }
      fprintf( out, "\n");
    }
  }
  else if ( format == "Full" ) {
   for ( int i=0; i<_origSize; i++ ) {
     for ( int j=0; j<_destSize; j++ ) {
       double val = getEntry( i, j );
       fprintf( out, " %10.8f", val );
     }
     fprintf( out, "\n");
    }
  }
  else if ( format == "R" ) {

    fprintf( out, "%s", toString("R").c_str() );
    fprintf( out, "\n" );

  }
  else if ( format == "SCILAB" ) {
    fprintf( out, "%s", toString("SCILAB").c_str() );

    // Scilab full format. Not complete.
 /*   fprintf( out, "full( [ ");
    double val;
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_destSize-1; j++ ) {
        val = getEntry( i, j );
        fprintf( out, " %10.8f,", val );
      }
      if ( i < _origSize-1 ) {
        val = getEntry( i, _origSize-1 );
        fprintf( out, " %10.8f;", val );
      }
      else {
        fprintf( out, " %10.8f]", val );
      }
    }
*/
  }
  else if ( format == "Matlab" ) {
    // Matlab sparse format
    // First print array of row indices
    fprintf( out, "theRows = [");
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_nbElts[i]; j++ ) {
        fprintf( out, " %4d", i );
      }
    }
    fprintf( out, "]';\n" );
    // Next print array of column indices
    fprintf( out, "theColumns = [");
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_nbElts[i]; j++ ) {
        fprintf( out, " %4d", _elts[i][j] );
      }
    }
    fprintf( out, "]';\n" );
    // Finally print array of values
    fprintf( out, "theValues = [");
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_nbElts[i]; j++ ) {
        fprintf( out, " %12f", _vals[i][j] );
      }
    }
    fprintf( out, "]';\n" );
    // To be used with "sparse( theRows, theColumns, theValues );
  }
  else {
    fprintf( stderr, "Warning in sparseMatrix::write(): format %s not recognized. Ignored.",
             format.c_str() );
  }

  free( destinations );
  free( values );

}

// method to consolidate transitions from some state
int sparseMatrix::consolidate(int i, int* destinations, double* values)
{
  int nbTrans = 0;

  for ( int e = 0; e < _nbElts[i]; e++ ) values[e] = 0.0;

  // straighforward insertion. No sorting, no shifting.
  for ( int e = 0; e < _nbElts[i]; e++ ) {
    int dest = _elts[i][e];
    bool isNew = true;
    int location = nbTrans;
    for ( int j = 0; isNew && ( j < nbTrans ); j++ ) {
      isNew = ( dest != destinations[j] );
      location = j;
    }
    if ( isNew ) {
      location = nbTrans;
      nbTrans++;
      destinations[location] = dest;
      values[location] = _vals[i][e];
    }
    else {
      values[location] += _vals[i][e];
    }
  }

  return nbTrans;
}

// method to normalize (sort rows) the matrix
void sparseMatrix::normalize()
{
  // Arrays necessary for consolidating transition probabilities.
  int maxNbElts = 0;
  for ( int i = 0; i < _origSize; i++ )
    if ( _nbElts[i] > maxNbElts ) maxNbElts = _nbElts[i];
  int* destinations = (int*)malloc( maxNbElts*sizeof(int) );
  double* values = (double*)malloc( maxNbElts*sizeof(double) );

  for ( int i = 0; i < _origSize; i++ ) {
    // normalize row i: consolidate then sort
    int nbTrans = consolidate( i, destinations, values );
    // bubble sort
    for ( int j = 0; j < nbTrans; j++ ) {
      for ( int k = j; k < nbTrans-1; k++ ) {
        if ( destinations[k+1] < destinations[k] ) {
          int d = destinations[k+1];
          destinations[k+1] = destinations[k];
          destinations[k] = d;
          double v = values[k+1];
          values[k+1] = values[k];
          values[k] = v;
        }
      }
    }
    // replace structures
    free( _elts[i] );
    free( _vals[i] );
    _nbElts[i] = nbTrans;
    _elts[i] = (int*) malloc( _nbElts[i] * sizeof(int) );
    _vals[i] = (double*) malloc( _nbElts[i] * sizeof(double) );
    for ( int j = 0; j < _nbElts[i]; j++ ) {
      _elts[i][j] = destinations[j];
      _vals[i][j] = values[j];
    }
  }

  free( destinations );
  free( values );
}
