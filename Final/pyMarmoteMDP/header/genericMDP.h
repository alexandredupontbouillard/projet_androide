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

#ifndef genericMDP_H
#define genericMDP_H

#include <string>
#include <vector>
#include "solutionMDP.h"
#include "Set/marmoteSet.h"
#include "transitionStructure/sparseMatrix.h"


/**
 * @brief Class genericMDP header : implementation of an abstract MDP class that represents an MDP object and can be inherited by other classes.
 * @author Hyon, lip6
 * @date 18 feb 2019
 * @version 1.4
 *
 * This class is generic, it contains the general elements of MDP: MDP type : time type, optimization rule, state space, action space,
 * action space, list of transition probabilities, list of transition rewards.
 *
 */
class genericMDP
{
public:
        /**
   	* @brief virtual function providing interface framework (printing function to stdout).
        */
	virtual void writeMDP();

	/**
   	* @brief Constructor to create genericMDP object.
   	* @author Hyon
	* @version 3
	* @date jan 2018
   	* @param t string : type_t value : "discrete" , "continuous".
	* @param r  string : type_r rule value can be  "min" , "max".
	* @param states marmoteSet : State space
	* @param actions marmoteSet : Action space
	* @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
	* @param rews sparseMatrix : reward structure : sparseMatrix object (state, action) an entry of the matrix is the reward for state line and  action column.
    * @return none.
    *
    * trans : transition structures vector of sparseMatrix objects an entry of the vector is a transition structure matrix from state to state
    * for the action at vector index.
    *
    * reward structure :  (state, action) an entry of the sparseMatrix is the reward for state line and  action column.
   	*/
	genericMDP(std::string t, std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, sparseMatrix* rews);

	/**
   	* @brief Constructor to create genericMDP object with a vector of transition structure.
   	* @author Hyon
	* @version 3
	* @date jan 2018
   	* @param t : string, type_t values : "discrete" , "continuous".
	* @param r : string, type_r rule can be  "min" , "max".
	* @param states : marmoteSet. State space
	* @param actions : marmoteSet. Action space
	* @param trans : vector<sparseMatrix*> transition structures vector of sparseMatrix.
	* @param rews : vector<sparseMatrix*> : vector of sparse matrix objects an entry of the vector is the reward structure from state to state for the action.
    * @return none.
    *
    * this second constructor take an other form  reward description.
    * It corresponds with the theoretical model where reward depends on transition r(i,a,j).
    * Then cost_perStage method is used to be consistent with the reward attribute
    *
    * trans : transition structures vector of sparseMatrix objects an entry of the vector is a transition structure matrix from state to state
    * for the action at vector index.
    *
    * rews : vector of sparse matrix objects an entry of the vector is the reward structure from state to state when trigger action (action is the index of the vector).
*/

	genericMDP(std::string t, std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, std::vector<sparseMatrix*> rews);

	//Destructor.
	virtual ~genericMDP();

    /**
   	* @brief A function to calculate costs (rewards) per stage.
   	* @author Abood Mourad.
    * @version 2
    * @date jan 2016
    * @param rews vector<sparseMatrix*> : for reward vector of transitionStructure objects (state, action, state), one reward matrix for each action.
    * @return sparseMatrix pointer.
    *
    * this method  transforms this vector in a transition matrix where first entry is the initial state and second entry is the action.
   */
    sparseMatrix* cost_perStage(std::vector<sparseMatrix*> rews);

    /**
      *@brief a function to delete reward when it has been created during the building of the object
      *@author EH
      *@date mar 2018
    */
    void clearRew();

    /**
    * @brief A function to solve (discrete time) MDP using value iteration algorithm.
    * @author Hyon
    * @version 3
    * @date jan 2018
    * @param epsilon double precision of the solution
    * @param maxIter int : the maximum number of iterations.
    * @return solutionMDP object.
   */
   virtual solutionMDP* valueIteration(double epsilon, int maxIter)=0;
   
   
   /**
    * @brief A function to solve (discrete time) MDP using value iteration algorithm with Gauss Seidel improvement.
    * @author Hyon
    * @version 1
    * @date feb 2019
    * @param epsilon double precision of the solution
    * @param maxIter int : the maximum number of iterations.
    * @return solutionMDP object.
   */
   virtual solutionMDP* valueIterationGS(double epsilon, int maxIter)=0;

    /**
    * @brief A function to solve (discrete time) MDP using modified policy iteration algorithm.
    * @author Hyon
    * @version 3
    * @date jan 2018
    * @param epsilon double.
    * @param maxIter int : the maximum number of iterations.
    * @param delta precision in the inner loop
    * @param moyIter nb iter max in the inner loop
    * @return solutionMDP object.
    *
    *
    */
	virtual solutionMDP* policyIterationModified(double epsilon, int maxIter, double delta, int moyIter)=0;


	/**
   	* @brief A function to solve (discrete time) MDP using policy iteration algorithm.
   	* @author Hyon
	* @version 3
	* @date jan 2018
	* @param maxIter int : the maximum number of iterations.
   	* @return solutionMDP object.
   	*/
	virtual solutionMDP* policyIteration(int maxIter) =0;

  /**
  * @brief A function to evaluate the cost of a policy with iteration of the power.
  * @author Hyon
  * @version 3
  * @date jan 2018
  * @param policy solutionMDP* : pointer to a solution that contains the optimal policy.
  * @param maxIter int : the maximum number of iterations.
  * @param epsilon double.
  * @return double.
  */
  virtual double* policyCost(solutionMDP *policy,double epsilon,int maxIter)=0;


protected:
	//Specific variables of the class:

	std::string type_t; 					/**< MDP type: "discrete" , "continuous". */
	std::string type_r;					/**< MDP rule: "min" , "max". */
	marmoteSet *stateSpace;					/**< State space represented by marmoteSet object. */
	marmoteSet *actionSpace;				/**< Action space represented by marmoteSet object. */
	std::vector<sparseMatrix*> transitions;       		/**< List (vector) of transitionStructure (one transition matrix for each action). */
	sparseMatrix* rewards;			         	/**< a rewards (costs) structure (row : state, column : action) */
};


#endif // genericMDP_H
