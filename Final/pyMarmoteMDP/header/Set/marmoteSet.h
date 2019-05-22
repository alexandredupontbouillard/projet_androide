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

/* marmoteSet.h ---
 *
 * Author: Alain Jean-Marie
 */

/* Commentary:
 *
 */

/* Change log:
 *
 */


#ifndef MARMOTESET_H
#define MARMOTESET_H

#include "../marmoteConstants.h"
#include <stdio.h>

/**
 * @brief The mother class representing abstract sets.
 *
 */
class marmoteSet {

 public:
  /**
   * @brief The different composition types for sets
   * - SIMPLE: the "atomic" typedef
   * - PRODUCT: a cartesian product of sets
   * - UNION: a disjoint union of subsets
   */
  enum opType { UNION, PRODUCT, SIMPLE };

 protected:
  opType _type; /**< type of the composition */

  int _nbDimensions; /**< number of dimensions if this is a product */
  int _nbZones; /**< number of subsets if this is a union */
  long int _cardinal; /**< cardinal of the set */
  marmoteSet** _zone; /**< array of zones, in case of a union */
  marmoteSet** _dimension; /**< array of dimensions, in case of a product */
  int* _stateBuffer; /**< a buffer to store the representation of a state as a vector if integers */

protected:
  int* _dimOffset; /**< table used for products to keep track of nested dimensions */
  int* _idxOffset; /**< table used for products to compute index offsets in dimensions */
  int _totNbDims; /**< total number of dimensions of the vector representation */
  int* _zeroState; /**< the representation of the initial state (static) */

public:
/**
 * @brief
 *
 */
  marmoteSet();
/**
 * @brief Constructor for composite sets, from smaller sets.
 *
 * @param list the list of elements in the composite
 * @param nb the number of elements in the composite
 * @param t the type of composite: UNION or PRODUCT
 */
  marmoteSet( marmoteSet **list, int nb, opType t );
  /**
   * Destructor for general sets. Essentially for composite structures.
   */
  virtual ~marmoteSet();

public:
  /**
   * @brief Cardinal of the set.
   *
   * @return an integer, the cardinal of the set, or INFINITE_STATE_SPACE_SIZE
   */
  virtual long int cardinal();
  /**
   * @brief Test if the set is finite.
   *
   * @return true if the set is finite, false if not.
   */
  virtual bool isFinite();
  /**
   * @brief Test if the set is a simple, elementary set.
   *
   * @return true if the set is elementary, false otherwise.
   */
  bool isSimple() { return _type == SIMPLE; };
  /**
   * @brief Test if the set is a composite, product set.
   *
   * @return true if the set is a product, false otherwise.
   */
  bool isUnion() { return _type == UNION; };
  /**
   * @brief Test if the set is a composite, product set.
   *
   * @return true if the set is a product, false otherwise.
   */
  bool isProduct() { return _type == PRODUCT; };
  /**
   * @brief Read accessor to the total number of dimensions.
   * @return The total number of dimensions of a vector representation.
   */
  int totNbDims() { return _totNbDims; };

 public:
  /**
   * @brief Enumerates the set: lists all elements of the set.
   * PROBLEM: enumerate where? in a file? to be specified more precisely.
   */
  virtual void enumerate();

  /**
   * @brief Tests the index() and decodeState() methods by enumerating all
   * states. Results written to standard output.
   */
  const void test_index_decode();

  /**
   * @brief Procedure that initializes the state buffer to the first state in the set
   * @param buffer the buffer to be initialized. Must have been allocated before.
   */
  virtual void firstState(int* buffer);
  /**
   * @brief Procedure to compute the state following a given state in the set order.
   * The new state is computed in place: the original state is erased. By convention,
   * the state following the last state is the first state.
   * @param buffer the reference state
   */
  virtual void nextState(int* buffer);
  /**
   * @brief Procedure that converts an index into a state. The state is written in
   * the buffer provided. The array must have been allocated before.
   * @param index index of the state to be decoded
   * @param buffer state buffer containing the resulting state
   */
  virtual void decodeState(int index, int* buffer);

  /**
   * @brief Function that computes the number (index) of some state in the order of the set.
   * The base class provides an implementation that uses solely firstState(), nextState(),
   * and state comparison. For efficiency, it is advised to re-implement this
   * method for derived classes.
   * @param buffer state buffer containing the state
   * @return the index of the state
   */
  virtual int index(int* buffer);

  /**
   * @brief Function that tests whether a state is the first state (state zero) or not.
   * @param buffer the state to be tested
   * @return true if the state is zero, false otherwise.
   */
  virtual bool isZero(int* buffer);
  /**
   * @brief Procedure to print a state, given by its full description, to some file descriptor.
   * @param out the file descriptor to be used
   * @param buffer the state to be printed
   */
  virtual void printState(FILE* out, int* buffer);
  /**
   * @brief Procedure to print a state, given by its index, to some file descriptor.
   * @param out the file descriptor to be used
   * @param index index of the state to be printed
   */
  virtual void printState(FILE* out, int index);
};

#endif
