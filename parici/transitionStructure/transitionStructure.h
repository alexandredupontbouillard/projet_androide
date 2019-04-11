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

#ifndef TRANSITIONSTRUCTURE_H
#define TRANSITIONSTRUCTURE_H

#include "../marmoteConstants.h"
#include "../Distribution/discreteDistribution.h"
#include <stdio.h>
#include <string>

/**
 * @brief Abstract class for transition structures. These are structures which describe
 * transitions to one state to another one, to which is attached a numeric label. Typical
 * instances should be one-step transition matrices of discrete-time Markov chains, and
 * infinitesimal generators of continuous-time Markov chains.
 * It is also possible that the origin state space and the destination state space are different.
 * @author Alain Jean-Marie and Issam Rabhi
 *
 */
class transitionStructure
{
 protected:
  timeType _type; /**< @brief the time type of the structure: discrete or continuous. */
                  /**< This determines the nature of the numeric label: probabilities or rates. */
  long int _origSize; /**< @brief size of the origin state space. */
  long int _destSize; /**< @brief size of the destination state space. */
  double _uniformizationRate; /**< @brief value related to the uniformization procedure. */

 public:
  // constructors
  // transitionStructure(int);
  /**
   * @brief Destructor of the class.
   *
   */
  virtual ~transitionStructure() {};

 public:
  //accessors
  /**
   * @brief Read accessor for the size of the state space. This is the
   * origin state space by default.
   *
   * @return the size of the origin state space
   */
  int size() { return _origSize; };
  /**
   * @brief Read accessor for the size of the origin state space.
   *
   * @return the size of the origin state space
   */
  int origSize() { return _origSize; };
  /**
   * @brief Read accessor for the size of the destination state space.
   *
   * @return the size of the destination state space
   */
  int destSize() { return _destSize; };
  /**
   * @brief Read accessor for the time type
   *
   * @return the time type
   */
  timeType type() { return _type; };
  /**
   * @brief Read accessor for the uniformization rate. Relevant mostly for discrete-time
   * structures created by uniformization of a continuous-time one.
   *
   * @return the uniformization rate
   */
  double uniformizationRate() { return _uniformizationRate; };
  /**
   * @brief Write accessor for the time type
   *
   * @param t the time type to be set
   */
  void setType(timeType t) { _type = t; };
  /**
   * @brief Write accessor for the uniformization rate.
   *
   * @param rate the rate to be set
   */
  void setUniformizationRate(double rate) { _uniformizationRate = rate; };
  /**
   * @brief Method to set the value associated with some transition.
   *
   * @param i the origin state
   * @param j the destination state
   * @param val the value attached to the transition
   * @return true if the operation was successful, false otherwise (out of bounds; wrong numeric value)
   */
  virtual bool setEntry(int i, int j, double val) = 0;
  /**
   * @brief Method to get the value associated with some transition. When state parameters
   * are out of bounds, the returned value should be 0.
   *
   * @param i the origin state
   * @param j the destination state
   * @return the value attached to the transition (i,j)
   */
  virtual double getEntry(int i, int j) = 0;
  /**
   * @brief Method to get the number of non-zero entries in a transition from some state.
   * This can be seen as the number of actually possible transitions from that state.
   *
   * @param i the origin state
   * @return int the number of possible transitions
   */
  virtual int getNbElts(int i) = 0;
  /**
   * @brief Method to get the @b number of the state corresponding to transition number k
   * in the list of possible transitions from some state i.
   *
   * @param i the origin state
   * @param k the index of transition from state i
   * @return the destination state corresponding to the k-th possible transition from state i
   */
  virtual int getCol(int i, int k) = 0;
  /**
   * @brief Method to get the @b value attached to transition number k
   * in the list of possible transitions from some state i.
   *
   * @param i the origin state
   * @param k the index of the transition from state i
   * @return the value attached to the k-th possible transition from state i
   */
  virtual double getEntryByCol(int i, int k) = 0;
  /**
   * @brief Method to get the transition from some state as a probability distribution.
   *
   * @param i the origin state
   * @return a discrete distribution object
   */
  virtual discreteDistribution* getTransDistrib(int i) = 0; // transitions from some state and their probas

  /**
   * @brief Reading from a file, and adding an element to the matrix.
   * The field must be in the form "row column value" with blank spaces
   * as separators
   * @author Alain Jean-Marie
   * @param input is the input flow from which the data is read
   * @return true if successful, false otherwise
   */
  bool readEntry(FILE *input);
  /**
   * @brief Calculating the sum of values corresponding to transitions from some state
   *
   * @param i the state
   * @return the sum
   */
  virtual double rowSum(int i) = 0;

 public:
  /**
   * @brief Copying a transition structure.
   *
   * @return transitionStructure
   */
  virtual transitionStructure* copy() = 0;
  /**
   * @brief Uniformizing a transition structure. The structure should be of continuous time type.
   * The resulting one will be of discrete time type. If uniformization fails, a NULL pointer
   * should be returned. If the origin structure is already of discrete-time type, a copy
   * should be returned.
   *
   * @return a discrete-time transition structure
   */
  virtual transitionStructure* uniformize() = 0;
  /**
   * @brief Embedding a discrete-time transition structure at jump times of a
   * transition structure. The structure should be of continuous time type.
   * The resulting one will be of discrete time type. If embedding fails, a NULL pointer
   * should be returned. If the origin structure is already of discrete-time type, a copy
   * should be returned.
   *
   * @return a discrete-time transition structure
   */
  virtual transitionStructure* embed() = 0;
  /**
   * @brief Computing the action of the transition structure on some measure, the measure
   * being represented as a vector of real numbers.
   * This corresponds to the multiplication vector/matrix, the row vector being interpreted
   * as a signed measure.
   *
   * @param d the measure to evaluate
   * @param res the resulting measure
   */
  virtual void evaluateMeasure(double* d,double* res) {
    fprintf( stderr, "Warning in transitionStructure::evaluateDistribution: not implemented.\n" ); }
  /**
   * @brief Computing the action of the transition structure on some probability distribution.
   * This version with discreteDistribution objects uses evaluateMeasure(double*,double*).
   *
   * @param d the distribution to evaluate
   * @return the distribution resulting from the evaluation
   */
  discreteDistribution* evaluateMeasure(discreteDistribution* d);
  /**
   * @brief Computing the action of the transition structure on some vector of values.
   * This corresponds to the multiplication matrix/vector, the column vector being interpreted
   * as a vector of values attached to the states.
   *
   * @param v the vector of values to evaluate
   * @param res the resulting vector
   */
  virtual void evaluateValue(double* v, double* res) = 0;

protected:
  // specific technical methods
  /**
   * @brief Method to consolidate (aggregate) transition probabilities from a given state.
   * The method returns the number of different destination states, and modifies the arrays
   * which will containt the lists of destinations and the corresponding probabilities.
   * These arrays must have been allocated beforehand, with a size large enough to handle
   * safely all possibilities.
   * @param i the origin state
   * @param destinations the list of destinations
   * @param values the list of probabilities
   * @return the number of different destinations
   */
  int consolidate(int i, int* destinations, double* values);

public:
  /**
   * @brief File output method for a transition structure
   *
   * @param out the file descriptor to which the structure should be written.
   * @param format the format/language to be used.
   */
  void write(FILE* out, std::string format);

  /**
   * @brief String serialization method for a transition structure. Supported formats are
   * XBORNE (Rii variant: by increasing row and increasing columns), MARCA, Matrix-Market sparse
   * and full, Ers, Maple, R, SCILAB, Full, and Matlab.
   * @author Issam Rabhi (R language), Alain Jean-Marie (XBORNE language)
   * @param format the format/language to be used.
   */
  std::string toString(std::string format);

};

#endif // TRANSITIONSTRUCTURE_H
