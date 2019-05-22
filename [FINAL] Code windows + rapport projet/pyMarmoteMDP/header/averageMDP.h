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

Copyright 2019 EMmanuel Hyon, Alain Jean-Marie*/

#ifndef averageMDP_H
#define averageMDP_H


#include "genericMDP.h"
#include "feedbackSolutionMDP.h"


/**
 * @brief Class averageMDP header: definition of an infinite horizon average MDP class.
 * @author Emmanuel Hyon.
 * @version 1.1
 * @date january 2019
 *
 * This class is inherited from the abstract class genericMDP.
 */

class averageMDP : public genericMDP
{
public:

    /**
    * @brief function providing printing function to stdout).
    */
    void writeMDP();

    /**
       * @brief Constructor to create averageMDP object.
       * @author EH
       * @version 0.91
       * @date jan 2019
       * @param r  string : type_r rule value can be  "min" , "max".
       * @param states marmoteSet : State space
       * @param actions marmoteSet : Action space
       * @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
       * @param rews sparseMatrix : reward structure : sparseMatrix object (state, action) an * entry of the matrix is the reward for state indexed by line and  action indexed by  * column.
       * @return none.
       *
       * trans : transition structures vector of sparseMatrix objects an entry of the vector * is a transition structure matrix from state to state
       * for the action at vector index.
       *
       * reward structure :  (state, action) an entry of the sparseMatrix is the reward for   * state in line and action in column.

    */
    averageMDP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, sparseMatrix* rews);

    /**
       * @brief Constructor to create averageMDP object.
       * @author EH
       * @version 0.9.1
       * @date jan 2019
       * @param r  string : type_r rule value can be  "min" , "max".
       * @param states marmoteSet : State space
       * @param actions marmoteSet : Action space
       * @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
       * @param rews vector<sparseMatrix*> vector of sparse matrix objects an entry of the
       * vector is the reward structure from state to state for the action.
       * @return none.
       *
       * This second constructor take an other form of reward description. It corresponds
       * with the theoretical model where reward depends on transition and reached state
       * r(i,a,j)
       * Then cost_perStage method is used to be consistent with the reward attribute
       *
       * trans : transition structures vector of sparseMatrix objects an entry of the vector
       * is a transition structure matrix from state to state
       * for the action at vector index.
       *
       * rews vector of sparse matrix objects: an entry of the vector is the reward structure
       * from state to state when trigger action (action is the index of the vector).
       *
    */
    averageMDP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, std::vector<sparseMatrix*> rews);

    /**
    * @brief destructor to delete averageMDP object.
    * @author EH
    */
    ~averageMDP(){};



  /**
  * @brief A function to change the index of the state used in value iteration
  * @author Hyon
  * @version 3
  * @date juil 2018
  * @param index int the index of the new
  */
  void changeIndex(int index){indexEtatSpecifique=index;}


  /**
  * @brief A function to solve (discrete time) MDP using value iteration algorithm.
  * @author Hyon
  * @version 0.1
  * @date jan 2019
  * @param epsilon double precision of the solution
  * @param maxIter int : the maximum number of iterations.
  * @return solutionMDP object.
  *
  * This is the value iteration algorithm which compute a Bellman Equation defined in part 8.5.1 of Puterman.
  *
  * In the value part of the solution. The average gain is placed. For an unichain problem   * the same value is given in the solution for any state.
  *
  */
  solutionMDP* valueIteration(double epsilon, int maxIter);
  
    /**
    * @brief A function to solve (discrete time) MDP using value iteration algorithm with Gauss Seidel improvement.
    * @author Hyon
    * @version 0
    * @date feb 2019
    * @param epsilon double precision of the solution
    * @param maxIter int : the maximum number of iterations.
    * @return solutionMDP object.
    * @warning not fully implemented
   */
   solutionMDP* valueIterationGS(double epsilon, int maxIter);

  
  
  /**
  * @brief A function to solve (discrete time) MDP using value relative iteration algorithm.
  * @author Hyon
  * @version 0.1
  * @date jan 2019
  * @param epsilon double precision of the solution
  * @param maxIter int : the maximum number of iterations.
  * @return solutionMDP object.
  *
  * This is the relative value iteration algorithm which compute a Bellman Equation including  * bias an average gain.
  *
  * In the value part of the solution. The average gain is placed. For an unichain problem   * the same value is given in the solution for any state.
  *
  */
  solutionMDP* relativeValueIteration(double epsilon, int maxIter);
  
  

  /**
    * @brief A function to solve (discrete time) MDP using policy iteration algorithm.
    * @author Hyon
    * @version 3
    * @date nov 2018
    * @param maxIter int : the maximum number of iterations.
    * @return solutionMDP object.
    * @warning not implemented
  */
  solutionMDP* policyIteration(int maxIter) ;

  /**
    * @brief A function to evaluate the average cost of a policy with iteration of the power.
    * @author Hyon
    * @version 0.1
    * @date Feb 2019
    * @param policy solutionMDP : object to get of the action
    * @param maxIter int : the maximum number of iterations.
    * @param epsilon double.
    * @return a pointer on a double which is the average cost in each state (identical in unichain model).
  */
  double* policyCost(solutionMDP *policy,double epsilon,int maxIter);

  /**
   * @brief A function to solve (discrete time) MDP using modified policy iteration algorithm.
   * @author Hyon
   * @version 0.1
   * @date Feb 2019
   * @param epsilon double the precision in the outer loop
   * @param maxIter int : the maximum number of iterations.
   * @param delta precision in the inner loop
   * @param maxInIter nb iter max in the inner loop
   * @return solutionMDP object.
   *
   * also called modified value iteration in Puterman's book
   *
   * In the value part of the solution. The average gain is placed.
   * For an unichain problem the same value is given in the solution for any state.
   *
  */
  solutionMDP* policyIterationModified(double epsilon, int maxIter, double delta, int maxInIter);


protected:
    //Specific variables of the class:
    double rho;                               /**< avergae value */
    std::string type_c;                       /**< MDP criteria: here "average" */
    int indexEtatSpecifique;                  /**< index de l'etat specifique pour calculer le cout moyen */
};

#endif // averageMDP_H
