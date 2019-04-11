/* Marmote and MarmoteMDP are free softwares: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Marmote is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Marmote. If not, see <http://www.gnu.org/licenses/>.

Copyright 2018 Emmanuel Hyon, Alain Jean-Marie, Abood Mourad
*/


#ifndef NONSTATIONARYSOLUTIONMDP_H
#define NONSTATIONARYSOLUTIONMDP_H

#include "solutionMDP.h"



/**
 * @brief Class nonStationarySolutionMDP: implementation of a nonStationarySolutionMDP class.
 * @author Abood Mourad, lip6 2016 Hyon lip6 2019.
 * @version 0.2
 * @date jan 2019
 *
 * This class is inherited from the class solutionMDP.
 * 
 */

class nonStationarySolutionMDP : public solutionMDP
{
public:

   /**
   * @brief Constructor to create feedbackSolution object.
   * @author EH
   * @version 0.2
   * @date jan 2019
   * @param h : int l'horizon du probleme fini
   * 
   * Constructeur qui initialise deux tableaux vides pour les deux attributs
   */
    nonStationarySolutionMDP(int h);

   /**
   * @brief destructor to delete feedbackSolution object.
   * @author AM EH
   * @version 0.2
   * @date jan 2019
   */
    ~nonStationarySolutionMDP();
    

    /**
    * @brief A function to print the solution object.
    * @author Abood Mourad Emmanuel Hyon.
    * @version 0.2
    * @date jan 2019 
    * 
    */
     void writeSolution();
     
   /**
   * @brief setter of the action vector at step t
   * @author EH
   * @date jan 2019
   * @version 0.1 
   * @param step : int the step for including the vector
   * @param a : int * . Pointeur sur un vecteur d'entiers de dimension le cardinal de l'espace d etat. 
   *   Chaque entree représente represente l'index de l'action à   
   *   effectuer quand on est dans l'etat dont l'index est la coordonnée. 
   * @return none.
   */
    void setAction(int step, int* a);

   /**
   * @brief setter of the value function vector at step t
   * @author EH
   * @date jan 2019
   * @version 0.1 
   * @param step : int the step for including the vector
   * @param t : double * Pointeur sur un vecteur de double de dimension le cardinal de l'espace d etat. 
   *  Chaque entree représente  la valeur de la fonction de valeur 
   *  dans l etat dont l'index est la coordonnée du vecteur
   * @return none.
   */
    void setValue(int step, double* t);

   /**
   * @brief setter of only one action at index at step t
   * @author EH
   * @date jan 2019
   * @version 0.1
   * @param step : int the step for including the vector
   * @param indice : int index de la valeur de l'action à modifier
   * @param value : int nouvelle valeur
   * @return none.
   */
    void setActionIndex(int step,int indice,int value);

   /**
   * @brief getter of only one action at index at step t
   * @author EH
   * @date jan 2019
   * @version 0.1
   * @param step : int the step for including the vector
   * @param indice : int index de la valeur de l'action à retourner
   * 
   * @return int the index of the action
   */
   int getActionIndex(int step, int indice);
     
     

//variables specifiques de la classe
protected:
    int horizon;       /**< the time horizon>*/
    int **action;      /**< An array whose columns are the index of actions representing the optimal policy. */
    double **value;    /**< An array whose columns are the value at each step. */
};

#endif // NONSTATIONARYSOLUTIONMDP_H




