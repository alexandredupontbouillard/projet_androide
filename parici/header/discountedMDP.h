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

Copyright 2018 Emmanuel Hyon, Alain Jean-Marie, Abood Mourad*/

#ifndef discountedMDP_H
#define discountedMDP_H


#include "genericMDP.h"
#include "feedbackSolutionMDP.h"


/**
 * @brief Class discountedMDP header: definition of an infinite horizon discounted MDP class.
 * @author Emmanuel Hyon
 * @version 3
 * @date dec 2018
 *
 * This class is inherited from the abstract class genericMDP.
 */

class discountedMDP : public genericMDP
{
public:
    /**
    * @brief function providing printing function to stdout.
    */
    void writeMDP();

    /**
      * @brief Constructor to create discountedMDP object.
      * @author EH
      * @version 1
      * @date dec 2018
      * @param r  string : type_r rule value can be  "min" , "max".
      * @param states marmoteSet : State space
      * @param actions marmoteSet : Action space
      * @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
      * @param rews sparseMatrix : reward structure : sparseMatrix object (state, action) an
      * entry of the matrix is the reward for state indexed by line and  action indexed by    * column.
      * @param b double : discount factor beta
      * @return none.
      *
      * trans : transition structures vector of sparseMatrix objects 
      * an entry of the vector is a transition structure matrix from state to state
      * for the action at vector index.
      *
      * reward structure :  (state, action) an entry of the sparseMatrix is the reward for
      * state in line and action in column.
    */
    discountedMDP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, sparseMatrix* rews, double b);

    /**
       * @brief Constructor to create discountedMDP object.
       * @author EH
       * @version 1
       * @date dec 2018
       * @param r  string : type_r rule value can be  "min" , "max".
       * @param states marmoteSet : State space
       * @param actions marmoteSet : Action space
       * @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
       * @param rews vector<sparseMatrix*> vector of sparse matrix objects an entry of the  
       * vector is the reward structure from state to state for the action.
       * @param b double : discount factor beta
       * @return none.
       *
       * This second constructor take an other form of reward description. It corresponds
       * with the theoretical model where reward depends on transition and reached state : 
       * r(i,a,j)
       * Then cost_perStage method is used to be consistent with the reward attribute
       *
       * trans : transition structures vector of sparseMatrix objects an entry of the vector * is a transition structure matrix from state to state
       * for the action at vector index.
       *
       * rews vector of sparse matrix objects an entry of the vector is the reward structure from state to state when trigger action (action is the index of the vector).
    */
    discountedMDP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, std::vector<sparseMatrix*> rews, double b);

    /**
    * @brief destructor to delete discountedMDP object.
    * @author EH
    */
    virtual ~discountedMDP();



  /**
  * @brief A function to solve (discrete time) MDP using value iteration algorithm.
  * @author Hyon
  * @version 3
  * @date jan 2018
  * @param epsilon double precision of the solution
  * @param maxIter int : the maximum number of iterations.
  * @return solutionMDP object.
  */
  solutionMDP* valueIteration(double epsilon, int maxIter);
  
  
   /**
    * @brief A function to solve (discrete time) MDP using value iteration algorithm with Gauss Seidel improvement.
    * @author Hyon
    * @version 1
    * @date feb 2019
    * @param epsilon double precision of the solution
    * @param maxIter int : the maximum number of iterations.
    * @return solutionMDP object.
    * @warning not fully implemented
   */
   solutionMDP* valueIterationGS(double epsilon, int maxIter);


  /**
   * @brief A function to solve (discrete time) MDP using modified policy iteration algorithm.
   * @author Hyon
   * @version 1.4
   * @date Feb 2019
   * @param epsilon double the precision in the outer loop
   * @param maxIter int : the maximum number of iterations.
   * @param delta precision in the inner loop
   * @param maxInIter nb iter max in the inner loop
   * @return solutionMDP object.
   *
   * also called hybrid value iteration in Powell
   *
  */
  solutionMDP* policyIterationModified(double epsilon, int maxIter, double delta, int maxInIter);

  /**
    * @brief A function to solve (discrete time) MDP using policy iteration algorithm.
    * @author Hyon
    * @version 3
    * @date jan 2018
    * @param maxIter int : the maximum number of iterations.
    * @return solutionMDP object.
    * @warning not fully implemented
  */
  solutionMDP* policyIteration(int maxIter) ;


  /**
    * @brief A function to evaluate the cost of a policy with iteration of the power.
    * @author Hyon
    * @version 1.1
    * @date april 2018
    * @param policy solutionMDP : object to get of the action
    * @param maxIter int : the maximum number of iterations.
    * @param epsilon double.
    * @return double.
  */
  double* policyCost(solutionMDP *policy,double epsilon,int maxIter);


  /**
    * @brief A function to evaluate the cost of a policy with iteration of the power and the scan of the set is done by index.
    * @author Hyon
    * @version 3
    * @date april 2018
    * @param policy solutionMDP : object to get of the action
    * @param maxIter int : the maximum number of iterations.
    * @param epsilon double.
    * @return double.
    *
    * For purpose of debug
  */
  double *policyCostbyIndex(solutionMDP *policy,double epsilon,int maxIter);

  protected:
  //Specific variables of the class:
  double beta;                  /**< discount factor */
  std::string type_c;           /**< MDP criteria: here "infinite discounted" */
};

#endif // discountedMDP_H
