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

#include "marmoteBox.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

// constructor for boxes with the SW corner at (0,...,0)
marmoteBox::marmoteBox(int nbDims, int *dimSize) : _nbDims(nbDims)
{
  // specific object variables
  _dimSize = (int*)malloc( _nbDims*sizeof(int) );
  _lVal = (int*)malloc( _nbDims*sizeof(int) );
  _uVal = (int*)malloc( _nbDims*sizeof(int) );

  _cardinal = 1;
  for ( int i = 0; i < _nbDims; i++ ) {
    _dimSize[i] = dimSize[i];
    _lVal[i] = 0;
    if ( _dimSize[i] == INFINITE_STATE_SPACE_SIZE ) {
       _uVal[i] = INFINITE_STATE_SPACE_SIZE;
    }
    else {
      _uVal[i] = _dimSize[i] - 1;
    }
    if ( ( _dimSize[i] == INFINITE_STATE_SPACE_SIZE )
         || ( _cardinal == INFINITE_STATE_SPACE_SIZE ) ) {
      _cardinal = INFINITE_STATE_SPACE_SIZE;
    }
    else {
      _cardinal *= _dimSize[i];
    }
  }

  // creating state buffer
  _stateBuffer = (int*)calloc( _nbDims, sizeof(int) );

  // initializing zero state to (0,0,...,0)
  _zeroState = (int*)calloc( _nbDims, sizeof(int) );

  // class variables
  _totNbDims = nbDims;

}

// constructor for general boxes
marmoteBox::marmoteBox(int nbDims, int *lower, int* upper)
{
  _nbDims = nbDims;
  _dimSize = (int*)malloc( _nbDims*sizeof(int) );
  _lVal = (int*)malloc( _nbDims*sizeof(int) );
  _uVal = (int*)malloc( _nbDims*sizeof(int) );

  _cardinal = 1;
  for ( int i = 0; i < _nbDims; i++ ) {
    _lVal[i] = lower[i];
    _uVal[i] = upper[i];
    if ( ( _lVal[i] != INFINITE_STATE_SPACE_SIZE ) && ( _uVal[i] < _lVal[i] ) ) {
      cerr << "Warning in marmoteBox(): lower value '" << _lVal[i]
              << "' larger than upper value '" << _uVal[i]
                 << "' in dimension #" << i << ". Interval ["
                 << _lVal[i] << ".." << _lVal[i]
                    << "] assumed." << endl;    _lVal[i] = 0;
      _uVal[i] = _lVal[i];
     }
    if ( _uVal[i] == INFINITE_STATE_SPACE_SIZE ) {
      _dimSize[i] = INFINITE_STATE_SPACE_SIZE;
    }
    else {
      _dimSize[i] = _uVal[i] - _lVal[i] + 1;
    }
    if ( ( _dimSize[i] == INFINITE_STATE_SPACE_SIZE )
         || ( _cardinal == INFINITE_STATE_SPACE_SIZE ) ) {
      _cardinal = INFINITE_STATE_SPACE_SIZE;
    }
    else {
      _cardinal *= _dimSize[i];
    }
  }

  // initializing zero state to lower values
  _zeroState = (int*)malloc( _nbDims*sizeof(int) );
  for ( int i = 0; i < _nbDims; i++ ) {
    _zeroState[i] = _lVal[i];
  }

  // class variables
  _totNbDims = nbDims;
}

// destructor
marmoteBox::~marmoteBox()
{
  free( _dimSize );
  free( _lVal );
  free( _uVal );
  free( _zeroState );
  free( _stateBuffer );
}

// finiteness test
bool marmoteBox::isFinite()
{
  bool res = true;

  for ( int i = 0; res && ( i < _nbDims); i++ ) {
    res = ( _dimSize[i] != INFINITE_STATE_SPACE_SIZE );
  }

  return res;
}

// test for the "zero" state
bool marmoteBox::isZero(int *buffer)
{
  bool res = true;

  for ( int d = 0; res && ( d < _nbDims); d++ ) {
    res = ( buffer[d] == _lVal[d] );
  }

  return res;
}

// initializing to the zero state
void marmoteBox::firstState(int *buffer)
{
  for ( int d = 0; d < _nbDims; d++ ) {
    buffer[d] = _lVal[d];
  }
}

// incrementing the state
void marmoteBox::nextState(int *buffer)
{
  buffer[_nbDims-1]++;
  for ( int d =_nbDims-1; ( d >= 0 ) && ( buffer[d] > _uVal[d] ); d-- ) {
    if ( d > 0 )
      buffer[d-1]++;
    buffer[d] = _lVal[d];
  }
}

// converting an index into an array
void marmoteBox::decodeState( int index, int* buf )
{
  int idx = index;
  for ( int d = _nbDims-1; d >= 0; d-- ) {
    buf[d] = _lVal[d] + idx % _dimSize[d];
    idx = idx / _dimSize[d];
  }
}

// index of some state
int marmoteBox::index(int *buf)
{
  int res = buf[0];
  for ( int d = 1; d < _nbDims; d++ ) {
    res = res * _dimSize[d] + buf[d];
  }

  return res;
}

// printing
void marmoteBox::printState(FILE *out, int *buffer)
{
  fprintf( out, "( " );
  for ( int d = 0; d < _nbDims-1; d++ ) {
    fprintf( out, "%3d, ", buffer[d] );
  };
  fprintf( out, "%3d", buffer[_nbDims-1] );
  fprintf( out, ")" );
}
