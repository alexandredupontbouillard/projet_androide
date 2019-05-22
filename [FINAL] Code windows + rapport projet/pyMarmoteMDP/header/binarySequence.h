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

#ifndef BINARYSEQUENCE_H
#define BINARYSEQUENCE_H

#include "marmoteSet.h"

/**
 * @brief The class representing binary sequences. These sets are cartesian products
 * of the basic set {0,1}. 
 *
 */

class binarySequence : public marmoteSet
{
private:
  int _dim; /**< number of elements in the sequence */

 public:
  /**
   * @brief Constructor for a binarySequence from the parameters n
   * @param n the number of elements
   */
  binarySequence(int n);
  /**
   * @brief Destructor.
   */
  ~binarySequence();
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

 private:
};

#endif // BINARYSEQUENCE_H
