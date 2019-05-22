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

#include "header/binarySimplex.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

// constructor for binary simplexes with parameters n and p
binarySimplex::binarySimplex(int n, int p)
{
  // specific object variables
  _nbPositions = n;
  _nbParticles = p;

  // creating state buffer
  _stateBuffer = (int*)calloc( _nbPositions, sizeof(int) );

  // initializing zero state to (1,1,...1,0,...,0)
  _zeroState = (int*)malloc( _nbPositions * sizeof(int) );
  for ( int i = 0; i < p; i++ )
    _zeroState[i] = 1;
  for ( int i = p; i < n; i++ )
    _zeroState[i] = 0;

  // class variables
  _totNbDims = n;

  // create binomial array
  _binom = (int**)malloc( (1+_nbPositions) * sizeof(int*) );
  for ( int i = 0; i <= _nbPositions; i++ )
    _binom[i] = (int*)calloc( 1+_nbPositions, sizeof(int) );
  
  _binom[0][0] = 1;
  for ( int j = 1; j <= _nbPositions; j++ )
    _binom[0][j] = 0;
  for ( int i = 1; i <= _nbPositions; i++ ) {
    _binom[i][0] = 1;
    for ( int j = 1; j <= _nbPositions; j++ ) {
      _binom[i][j] = _binom[i-1][j-1] + _binom[i-1][j];
    }
  }

  // compute cardinal as: n choose p
  //  _cardinal = 1L;
  // for ( int i = 1; i < _nbParticles; i++ ) {
  //   _cardinal = (_cardinal * ( n+1-i ) ) / i;
  // }
  _cardinal = _binom[_nbPositions][_nbParticles];
  
  // dump it
  // for ( int i = 0; i <= _nbPositions; i++ ) {
  //   for ( int j = 0; j <= _nbPositions; j++ ) {
  //     printf("%4d", _binom[i][j] );
  //   }
  //   printf("\n");
  // }
  
}

// destructor
binarySimplex::~binarySimplex()
{
  free( _zeroState );
  free( _stateBuffer );
  
  for ( int i = 0; i <= _nbPositions; i++ ) {
    free( _binom[i] );
  }
  free( _binom );
}

// finiteness test
bool binarySimplex::isFinite()
{
  bool res = true;

  return res;
}

// test for the "zero" state
bool binarySimplex::isZero(int *buffer)
{
  bool res = true;

  for ( int d = 0; res && ( d < _nbPositions ); d++ ) {
    res = ( buffer[d] == _zeroState[d] );
  }

  return res;
}

// initializing to the zero state
void binarySimplex::firstState(int *buffer)
{
  for ( int d = 0; d < _nbPositions; d++ ) {
    buffer[d] = _zeroState[d];
  }
}

// incrementing the state
void binarySimplex::nextState(int *buffer)
{
  // use standard implementation: increment index and decode
  int idx = binarySimplex::index(buffer);
  if ( idx < _cardinal-1 ) {
    binarySimplex::decodeState(idx+1,buffer);
  }
  else {
    binarySimplex::firstState(buffer);
  }
}

// converting an index into an array
void binarySimplex::decodeState( int index, int* buf )
{
  int idx = index;
  int part_left = _nbParticles;
  int v;
  
  for ( int d = 0; d < _nbPositions; d++ ) {
    if ( part_left == 0 ) {
      v = 0;
    }
    else {
      v = _binom[_nbPositions-1-d][part_left-1];
    }
    // printf("Decoding: position %2d, left particles = %d, idx = %d, seuil = %d\n",
    // d, part_left, idx, v );
    if ( idx < v ) {
      buf[d] = 1;
      part_left--;
    }
    else {
      buf[d] = 0;
      idx -= v;
    }
  }
}

// auxiliary recursive function to evaluate indices
int binarySimplex::idx(int n, int p, int* buf) {
  int res = 0;

  if ( p == 0 ) {
    res = 0;
  }
  else {
    if ( buf[0] == 0 ) {
      res = _binom[n-1][p-1] + binarySimplex::idx( n-1, p, &buf[1] );
    }
    else {
      res = binarySimplex::idx( n-1, p-1, &buf[1] );
    }
  }
  
  return res;
}

// index of some state
int binarySimplex::index(int *buf)
{
  int res = idx( _nbPositions, _nbParticles, buf );

  return res;
}

// printing
void binarySimplex::printState(FILE *out, int *buffer)
{
  fprintf( out, "(" );
  for ( int d = 0; d < _nbPositions; d++ ) {
    fprintf( out, "%2d", buffer[d] );
  };
  fprintf( out, " )" );
}
