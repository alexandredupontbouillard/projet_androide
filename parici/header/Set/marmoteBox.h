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

#ifndef MARMOTEBOX_H
#define MARMOTEBOX_H

#include "marmoteSet.h"

/**
 * @brief The class representing "rectangular" sets. Boxes are cartesian products
 * of intervals.
 *
 */

class marmoteBox : public marmoteSet
{
private:
  int _nbDims; /**< number of dimensions of the box */
  int* _dimSize; /**< number of states in each dimension */
  int* _lVal; /**< array of lower values of intervals in each dimension */
  int* _uVal; /**< array of upper values of intervals in each dimension */

public:
  /**
   * @brief Constructor for a marmoteBox from the arrays of sizes in each dimension.
   * By convention, the SW corner of the box is (0,...,0).
   * All arrays passed as parameters are copied.
   * @param nbDims the number of dimensions
   * @param dimSize array of size in each dimension (possibly INFINITE_STATE_SPACE_SIZE)
   */
  marmoteBox(int nbDims, int* dimSize);
  /**
   * @brief Constructor for general boxes, from the arrays of lower/upper values.
   * By convention, these lower/upper values belong to the box.
   * All arrays passed as parameters are copied.
   * @param nbDims the number of dimensions
   * @param lower array of lower values in each dimension
   * @param upper array of upper values in each dimension (possibly INFINITE_STATE_SPACE_SIZE)
   */
  marmoteBox(int nbDims, int *lower, int* upper);
  /**
   * @brief Destructor.
   */
  ~marmoteBox();
  /**
   * @brief Test whether the set is finite. Boxes are finite if and only if the
   * size in each dimension is finite.
   */
  bool isFinite();
  /**
   * @brief Tests that a state, given by its vector representation, is zero.
   * @param buffer the state to be tested
   * @return true if the state is 0, false otherwise.
   */
  bool isZero(int* buffer);
  /**
   * @brief Initializes some state buffer with the first state of the set.
   * @param buffer the buffer to be set.
   */
  void firstState(int* buffer);
  /**
   * @copydoc marmoteSet::nextState()
   */
  void nextState(int *buffer);
  /**
   * @brief utility to convert a state index into a state array
   * @author Alain Jean-Marie
   * @param index the state index
   * @param buf the state buffer to be filled
   */
  void decodeState(int index, int* buf);
  /**
   * @brief Utility to find the number of some state.
   * @param buf a state buffer
   * @return the index of the state buffer
   */
  int index(int* buf);
  /**
   * @brief Procedure for printing out a state.
   * @param out file descriptor of the stream to be used
   * @param buffer the state buffer to be printed
   */
  void printState(FILE* out, int *buffer);

};

#endif // MARMOTEBOX_H
