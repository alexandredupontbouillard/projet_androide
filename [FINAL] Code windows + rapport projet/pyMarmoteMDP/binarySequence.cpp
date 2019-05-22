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

Copyright 2016 Alain Jean-Marie, Jean-Michel Fourneau, Jean-Marc Vincent */

#include "header/binarySequence.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

// constructor for binary sequences with parameter n
binarySequence::binarySequence(int n): _dim(n)
{
  // creating state buffer
  _stateBuffer = (int*)calloc( _dim, sizeof(int) );

  // initializing zero state to (0,...,0)
  _zeroState = (int*)malloc( _dim * sizeof(int) );
  for ( int i = 0; i < _dim; i++ )
    _zeroState[i] = 0;

  // class variables
  _totNbDims = _dim;

  // compute cardinal as: 2**n
  _cardinal = 1;
  for ( int i = 1; i <= _dim; i++ ) _cardinal *= 2;
  
}

// destructor
binarySequence::~binarySequence()
{
  free( _zeroState );
  free( _stateBuffer );
}

// finiteness test
bool binarySequence::isFinite()
{
  bool res = true;

  return res;
}

// test for the "zero" state
bool binarySequence::isZero(int *buffer)
{
  bool res = true;

  for ( int d = 0; res && ( d < _dim ); d++ ) {
    res = ( buffer[d] == _zeroState[d] );
  }

  return res;
}

// initializing to the zero state
void binarySequence::firstState(int *buffer)
{
  for ( int d = 0; d < _dim; d++ ) {
    buffer[d] = _zeroState[d];
  }
}

// incrementing the state
void binarySequence::nextState(int *buffer)
{
  // use standard implementation: increment index and decode
  int idx = binarySequence::index(buffer);
  if ( idx < _cardinal-1 ) {
    binarySequence::decodeState(idx+1,buffer);
  }
  else {
    binarySequence::firstState(buffer);
  }
}

// converting an index into an array
void binarySequence::decodeState( int index, int* buf )
{
  int remainder = index;

  for ( int d = _dim-1; d >= 0; d-- ) {
    buf[d] = remainder % 2;
    remainder /= 2;
  }
}

// Index of some state. 
int binarySequence::index(int *buf)
{
  int res = 0;

  for ( int i = 0; i < _dim; i++ ) {
    res *= 2;
    res += buf[i];
  }

  return res;
}

// printing
void binarySequence::printState(FILE *out, int *buffer)
{
  fprintf( out, "(" );
  for ( int d = 0; d < _dim; d++ ) {
    fprintf( out, "%2d", buffer[d] );
  };
  fprintf( out, " )" );
}
