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

#ifndef EVENTMIXTURE_H
#define EVENTMIXTURE_H

#include "transitionStructure.h"
#include "sparseMatrix.h"
#include <string>

/**
 * @brief A class representing probabilistic mixtures of "events", that is, elementary 
 * transition structures.
 *
 */

class eventMixture : public transitionStructure
{
private:
  // specific attributes
  int _nbEvents; /**< number of different events */
  double* _eventProba; /**< probability of each event */
  std::string* _eventName; /**< names given to the events */
  int** _eventTransition; /**< table of state-to-state mappings corresponding to each transition */

public:
  // constructors, destructor
  /**
   * @brief Constructor from arrays.
   * @param size the number of states
   * @param nbEvents the number of different events
   * @param probas the array of probabilities. This array is <b>copied</b>
   * @param names the names of events. This array is <b>copied</b>
   * @param transitions the array of state-to-state mappings. This array is <b>not copied</b>.
   */
  eventMixture(int size, int nbEvents, double* probas, std::string* names, int** transitions );
  /**
   * @brief Contruction from a sparse matrix.
   * @param spMat the sparse matrix to convert into an event mixture
   */
  eventMixture(sparseMatrix* spMat);
  /**
   * @brief Destructor of the class.
   *
   */
  ~eventMixture();

public:
  // accessors
  /**
   * @brief Read accessor for the number of events
   * @return the number of events
   */
  int nbEvents() { return _nbEvents; }
  /**
   * @brief Read accessor for the probabilities associated to each event.
   * @param e the event index
   * @return the probability that the event is e
   */
  double eventProba(int e) { return _eventProba[e]; }

public:
  // implementation of virtual methods
  bool setEntry(int i, int j, double val);
  double getEntry(int i, int j);
  int getNbElts(int i);
  int getCol(int i, int k);
  double getEntryByCol(int i, int k);
  discreteDistribution* getTransDistrib(int i);
  double rowSum(int i);
  eventMixture* copy();
  eventMixture* uniformize();
  eventMixture* embed();
  /**
   * @copydoc transitionStructure::evaluateMeasure(double*,double*)
   */
  void evaluateMeasure(double* d, double* res) ;

  /**
   * @copydoc transitionStructure::evaluateMeasure(discreteDistribution*)
   */
  discreteDistribution* evaluateMeasure(discreteDistribution* d);

  /**
   * @brief Evaluation of the expected value after a transition, for all states.
   * Corresponds to matrix/vector multiplication. The result is stored in an array
   * that must be already allocated.
   * @author Alain Jean-Marie
   * @param v the value (column vector) to be multiplied
   * @param res the resulting vector
   */
  void evaluateValue(double* v, double* res);
   /**
   * @brief Evaluation of the expected value after a transition from some state.
   * Corresponds to matrix/vector multiplication for a specific state/line of the matrix.
   * @author Alain Jean-Marie
   * @param v the value (column vector) to be multiplied
   * @param stateIndex the index of the state for which to perform the multiplication
   * @return the value of the state
   * */
  double evaluateValueState(double* v, int stateIndex);
  /**
   * @brief Output method for the transition structure. Supported formats are:
   * XBORNE, MARCA, Ers, Maple and PSI3 (experimental)
   *
   * @param out the file descriptor to which the structure should be written.
   * @param format the format/language to be used
   */
  void write(FILE* out, std::string format);

private:
  // specific technical methods
  /**
   * @brief Method to consolidate (aggregate) transition probabilities from a given state.
   * The method returns the number of different destination states, and modifies the arrays
   * which will containt the lists of destinations and the corresponding probabilities.
   * These arrays must have been allocated beforehand, with a size at least the maximum number of
   * transitions from any state.
   * @param i the origin state
   * @param destinations the list of destinations
   * @param values the list of probabilities
   * @return the number of different destinations
   */
  int consolidate(int i, int* destinations, double* values);
};

#endif // EVENTMIXTURE_H
