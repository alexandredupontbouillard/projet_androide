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

#include "marmoteSet.h"
#include <stdlib.h>
#include <iostream>

// Standard constructor. Represents the empty set.
marmoteSet::marmoteSet()
{
  _type = SIMPLE;
  _cardinal = 0;
  _nbZones = 0;
  _nbDimensions = 0;
  _totNbDims = 0;
  _zone = (marmoteSet**)NULL;
  _dimension = (marmoteSet**)NULL;
  _zeroState = (int*)NULL;
  _stateBuffer = (int*)NULL;
  _dimOffset = (int*)NULL;
}

// constructor of composite sets from
marmoteSet::marmoteSet( marmoteSet **list, int nb, opType t )
{
  switch( t ) {
  case UNION:
    _type = UNION;
    _nbZones = nb;
    _totNbDims = 0;

    for ( int i = 0; i < _nbZones; i++ ) {
      _zone[i] = list[i];
      if ( _zone[i]->totNbDims() > _totNbDims ) {
        _totNbDims = _zone[i]->totNbDims();
      }
    }

    // initialize state buffer and "zero" state as zero state of the first
    // subset
    _stateBuffer = (int*)malloc( _totNbDims*sizeof(int) );
    _zeroState = (int*)malloc( _zone[0]->totNbDims()*sizeof(int) );
    _zone[0]->firstState( _zeroState );

    break;

  case PRODUCT:
    _type = PRODUCT;
    _nbDimensions = nb;
    _totNbDims = 0;
    
    _dimension = (marmoteSet**)malloc( _nbDimensions*sizeof(marmoteSet*) );
    _dimOffset = (int*)malloc( _nbDimensions*sizeof(int) );
    _dimOffset[0] = 0;
    for ( int i = 0; i < _nbDimensions; i++ ) {
      _dimension[i] = list[i];
      _totNbDims += _dimension[i]->totNbDims();
      if ( i < _nbDimensions-1 ) {
        _dimOffset[i+1] = _dimOffset[i] + _dimension[i]->totNbDims();
      }
    }

    // initialize state buffer and "zero" state with zero states in each
    // dimension
    _stateBuffer = (int*)malloc( _totNbDims*sizeof(int) );
    _zeroState = (int*)malloc( _totNbDims*sizeof(int) );
    for ( int i = 0; i < _nbDimensions; i++ ) {
      _dimension[i]->firstState( &(_zeroState[_dimOffset[i]]) );
    }

    break;

  case SIMPLE:
    fprintf( stderr, "Error: improper use of marmoteSet constructor with SIMPLE type.\n");
  }

  _cardinal = cardinal();

}

marmoteSet::~marmoteSet()
{
  // It is more complicated than expected to free internal structures.
  // Instantiations of marmoteSet of the "PRODUCT" or "UNION" type
  // do create their "_dimension" or "_zone" themselves and have to
  // free it themselves. I don't know how to make the difference.
  /*
  if ( _type == PRODUCT ) {
    free( _dimension );
  }
  else if ( _type == UNION ) {
    free( _zone );
  }
  else {
    // elementary structures free their private arrays themselves
  }
  */
}

long int marmoteSet::cardinal()
{
  long int res;

  if ( _type == PRODUCT ) {
    res = 1;
    for ( int i = 0; i < _nbDimensions; i++ ) {
      res = res * _dimension[i]->cardinal();
    }
  }
  else if ( _type == UNION ) {
    res = 0;
    for ( int i = 0; i < _nbZones; i++ ) {
      res = res + _zone[i]->cardinal();
    }
  }
  else {
    // must be atomic
    res = _cardinal;
  }

  return res;
}

bool marmoteSet::isFinite()
{
  bool res = true;

  if ( _type == PRODUCT ) {
    for ( int i = 0; i < _nbDimensions; i++ ) {
      res = res && _dimension[i]->isFinite();
    }
  }
  else if ( _type == UNION ) {
    for ( int i = 0; i < _nbZones; i++ ) {
      res = res && _zone[i]->isFinite();
    }
  }
  else {
    // must be atomic
    fprintf( stderr, "Error: in marmoteSet::isFinite(): elementary sets must have their own method. True assumed.\n" );
    
  }

  return res;
}

const void marmoteSet::test_index_decode()
{
  // use one's own buffer in order not to intefere with other functions
  int* testBuf = (int*)malloc( _totNbDims * sizeof(int) );
  
  for ( int i = 0; i < _cardinal; i++ ) {
    decodeState(i,testBuf);
    std::cout << "Idx " << i << " ";
    printState(stdout,testBuf);
    std::cout << " -> " << index(testBuf) << std::endl;
  }

  free(testBuf);
}

void marmoteSet::enumerate()
{

  switch( _type ) {
  case UNION:
    for ( int i = 0; i < _nbZones; i++ ) {
      _zone[i]->enumerate();
    }

  case PRODUCT:
    // lazy implementation for products...
  case SIMPLE:
    // must be atomic
    // selfEnumerate();
    firstState(_stateBuffer);
    do {
      printState( stdout, _stateBuffer );
      nextState( _stateBuffer );
    } while( !isZero( _stateBuffer ) );
  }

}

// test whether some state buffer is the zero state
bool marmoteSet::isZero(int* buffer)
{
  bool res = true;

  if ( _type == PRODUCT ) {
    int* ptr;
    for( int i = 0; res && ( i < _nbDimensions ); i++ ) {
      ptr = &(buffer[_dimOffset[i]]);
      res = _dimension[i]->isZero( ptr );
    }
  }
  else if ( _type == UNION ) {
    // no idea to which subset to ask...
  }
  else {
    // set is elementary
    fprintf( stderr, "Error: in marmoteSet::isZero(): elementary sets must have their own method. Ignored.\n" );
  }

  return res;
}

// set, in place, some state buffer to the zero state
void marmoteSet::firstState(int* buffer)
{
  switch ( _type ) {
  case marmoteSet::PRODUCT:
    int* ptr;
    for( int i = 0; i < _nbDimensions; i++ ) {
      ptr = &(buffer[_dimOffset[i]]);
      _dimension[i]->firstState( ptr );
    }
    break;

  case marmoteSet::UNION:
    _dimension[0]->firstState(buffer);
    break;

  case marmoteSet::SIMPLE:
    // set is elementary
    fprintf( stderr, "Error: in marmoteSet::firstState(): elementary sets must have their own method. Ignored.\n" );
    break;
  }
}

// increment, in place, the state buffer
void marmoteSet::nextState(int* buffer)
{
  switch ( _type ) {
  case marmoteSet::PRODUCT:
    // scan dimensions from right to left, and carry over when dimension
    // is back to 0
    int numDim;
    numDim = _nbDimensions-1;
    int* ptr;
    do {
      ptr = &(buffer[_dimOffset[numDim]]);
      _dimension[numDim]->nextState( ptr );
      numDim--;
    } while( ( numDim >= 0 ) && _dimension[numDim+1]->isZero( ptr ) );
    break;

  case marmoteSet::UNION:
    // in which subset is this state?? no idea!!
    break;

  case marmoteSet::SIMPLE:
    // set is elementary
    fprintf( stderr, "Error: in marmoteSet::nextState(): elementary sets must have their own method. Ignored.\n" );
  }
}

// converts an index into a state buffer
void marmoteSet::decodeState(int index, int* buffer)
{
  switch ( _type ) {
  case marmoteSet::PRODUCT:
    int idx;
    idx = index;
    
    for ( int d = _nbDimensions-1; d >= 0; d-- ) {
      int local_idx;
      local_idx = idx % _dimension[d]->cardinal();
      int* ptr;
      ptr = &(buffer[_dimOffset[d]]);
      _dimension[d]->decodeState(local_idx,ptr);

      idx = idx / _dimension[d]->cardinal();
    }
    break;
    
  case marmoteSet::UNION:
    // in which subset is this state?? no idea!!
    break;

  case marmoteSet::SIMPLE:
    // set is elementary
    fprintf( stderr, "Error: in marmoteSet::decodeState(): elementary sets must have their own method. Ignored.\n" );
  }
}

// converts a state into an index
int marmoteSet::index(int* buffer)
{
  int res = 0;
  
  switch( _type ) {
  case PRODUCT:
    int d;
    int* ptr;

    ptr = &(buffer[_dimOffset[0]]);
    res = _dimension[0]->index( ptr );
    for ( d = 1; d < _nbDimensions; d++ ) {
      ptr = &(buffer[_dimOffset[d]]);
      res = res * _dimension[d]->cardinal() + _dimension[d]->index( ptr );
    }

    break;

  default:
    // Just scans through the states until the right one is found.
    // Relies on a sound "nextState" method.
    firstState(_stateBuffer);
    bool found = false;
    for ( res = 0; !found && ( res < _cardinal); res++ ) {
      bool isEqual = true;
      for ( int j = 0; isEqual && ( j < _totNbDims ); j++ ) {
	isEqual = ( _stateBuffer[j] == buffer[j] );
      }
      found = isEqual;
      nextState(_stateBuffer);
    }
    
    if ( !found ) {
      res = -1;
    }
    else {
      res--;
    }
  }

  return res;
}

// printing a state as a vector
void marmoteSet::printState(FILE* out, int* buffer)
{
  switch ( _type ) {
  case marmoteSet::PRODUCT:
    fprintf( out, "(" );
    int* ptr;
    for ( int i = 0; i < _nbDimensions; i++ ) {
      ptr = &(buffer[_dimOffset[i]]);
      _dimension[i]->printState( out, ptr );
      if ( i != _nbDimensions-1 ) {
        fprintf( out, ";" );
      }
    }
    fprintf( out, ")" );
    break;

  case marmoteSet::UNION:
    // in which subset is this state?? no idea!!
    break;

  case marmoteSet::SIMPLE:
    // set is elementary
    fprintf( stderr, "Error: in marmoteSet::printState(): elementary sets must have their own method. Ignored.\n" );
  }
}

// printing a state as a vector, from its index number
void marmoteSet::printState(FILE *out, int index)
{
  decodeState( index, _stateBuffer );
  printState( out, _stateBuffer );
}
