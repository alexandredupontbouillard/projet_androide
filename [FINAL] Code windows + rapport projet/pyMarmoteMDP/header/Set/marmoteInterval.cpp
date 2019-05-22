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

#include "marmoteInterval.h"
#include <stdlib.h>

// constructor from interval ends
marmoteInterval::marmoteInterval( int min, int max )
{
  _min = min;
  _max = max;

  if ( _min <= _max ) {
    _cardinal = _max - _min + 1;
  }
  else {
    // empty interval
    _cardinal = 0;
  }

  // initialization of other mother class variables
  _type = SIMPLE;
  _nbZones = 0;
  _nbDimensions = 0;
  _totNbDims = 1;
  _zone = (marmoteSet**)NULL;
  _dimension = (marmoteSet**)NULL;
  _stateBuffer = (int*)malloc( 1 * sizeof(int) );
  _zeroState = (int*)malloc( 1 * sizeof(int) );
}

marmoteInterval::~marmoteInterval() {
  free( _stateBuffer );
  free( _zeroState );
}

// test for the "zero" state
bool marmoteInterval::isZero(int *buffer)
{
  bool res = (buffer[0] == _min );

  return res;
}

// initializing to the zero state
void marmoteInterval::firstState(int *buffer)
{
  buffer[0] = _min;
}

// incrementing the state
void marmoteInterval::nextState(int *buffer)
{
  buffer[0]++;
  if ( buffer[0] > _max )
    buffer[0] = _min;
}

// converting an index into an array
void marmoteInterval::decodeState(int index, int* buffer )
{
  buffer[0] = _min + index;
}

// index of some state
int marmoteInterval::index(int *buffer)
{
  int res = buffer[0] - _min;

  return res;
}

// printing
void marmoteInterval::printState(FILE *out, int *buffer)
{
  fprintf( out, " %4d", buffer[0] );
}

// enumeration
void marmoteInterval::enumerate()
{
  int buffer[1];

  for ( int i = _min; i <= _max; i++ ) {
    buffer[0] = i;
    printState( stdout, buffer );
  }
}
