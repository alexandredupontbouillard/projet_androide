
#ifndef discountedMDPLP_H
#define discountedMDPLP_H


#include "genericMDP.h"
#include "feedbackSolutionMDP.h"
#include "discountedMDP.h"

/**
 * @brief Class discountedMDP_LP header: definition of an infinite horizon discounted MDP class which would be solved by linear programming
 * @author Hyon, lip6.
 * @version 0.1
 * @date june 2018
 *
 * This class is inherited from the class discountedMDP. 
 */

class discountedMDP_LP : public discountedMDP
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
       * @date jan 2018
       * @param r  string : type_r rule value can be  "min" , "max".
       * @param states marmoteSet : State space
       * @param actions marmoteSet : Action space
       * @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
       * @param rews sparseMatrix : reward structure : sparseMatrix object (state, action) an entry of the matrix is the reward for action line and state column.
       * @param b double : discount factor beta
       * @return none.
       * 
       * trans : transition structures vector of sparseMatrix objects an entry of the vector is a transition structure matrix from state to state 
       * for the action at vector index.
       *
       * reward structure :  (state, action) an entry of the sparseMatrix is the reward for action line and state column.     

       */
       discountedMDP_LP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, sparseMatrix* rews, double b);

       /**
       * @brief Constructor to create discountedMDP object.
       * @author EH
       * @version 1
       * @date jan 2018
       * @param r  string : type_r rule value can be  "min" , "max".
       * @param states marmoteSet : State space
       * @param actions marmoteSet : Action space
       * @param trans vector<sparseMatrix*> : transition structures vector of sparseMatrix.
       * @param rews vector<sparseMatrix> vector of sparse matrix objects an entry of the vector is the reward structure from state to state for the action.
       * @param b double : discount factor beta
       * @return none.
       *
       * This second constructor take an other form  reward description. It corresponds with the theoretical model where reward depends on transition r(i,a,j) 
       * Then cost_perStage method is used to be consistent with the reward attribute
       *
       * trans : transition structures vector of sparseMatrix objects an entry of the vector is a transition structure matrix from state to state 
       * for the action at vector index.
       *
       * rews vector of sparse matrix objects an entry of the vector is the reward structure from state to state when trigger action (action is the index of the vector).

       */
       discountedMDP_LP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, std::vector<sparseMatrix*> rews, double b);  

       /**
       * @brief destructor to delete discountedMDP object.
       * @author EH
       */
       ~discountedMDP_LP();
       
       
       
       /**
       * @brief A function to solve (discrete time) MDP using value iteration algorithm.
       * @author Hyon
       * @version 3
       * @date jan 2018
       * @param epsilon double precision of the solution
       * @param maxIter int : the maximum number of iterations.
       * @return solutionMDP object.
       */
       solutionMDP* linearProgrammingSolve();

protected :   
    double ** constraints; /**< matrix for the constraints */
    double * objectives; /**< matrix for the objectives */
       
private:
    
    /**
       * @brief A function to build the linear object
       * @author Hyon
       * @version 0.1
       * @date june 2018
       * 
       * dedicated gurobi
       */
       void gurobiBuild(){};
    
};

#endif // discountedMDPLP_H
