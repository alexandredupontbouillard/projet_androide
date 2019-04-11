
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

Copyright 2019 Emmanuel Hyon, Alain Jean-Marie*/

/**
 * @brief Class totalRewardMDP source
 * @author Hyon, lip6
 * @date feb 2019
 * @version 1
 *
 * This .cpp file is for implementing functions of a total reward MDP class class
 */

#include <values.h>
#include <iostream>
#include "alglin.h"
#include "header/totalRewardMDP.h"

using namespace std;


// constructeurs

totalRewardMDP::totalRewardMDP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, sparseMatrix* rews)
: genericMDP("discrete",r,states,actions,trans,rews), type_c("total reward")
{
}


totalRewardMDP::totalRewardMDP(std::string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, vector<sparseMatrix*> rews)
: genericMDP(string("discrete"),r,states,actions,trans,rews), type_c("total reward")
{
}

// destructeur
totalRewardMDP::~totalRewardMDP(){
}

// other utilities function
void totalRewardMDP::writeMDP()
{
    genericMDP::writeMDP();
    printf("#############################################\n");
    printf("MDP Type : %s \n" , type_c.c_str());
}

//A function to solve infinite total reward  MDP using value iteration algorithm.
solutionMDP* totalRewardMDP::valueIteration(double epsilon, int maxIter)
{
  //Create necessary variables.
    int i;
    int a;
    int nbIter; // count the iteration number
    double optimal; //keep optimal value
    int argopt;
    int indiceEtat;
    int indiceAction;

    double distance;
    double gain;
    double res;

    double *newVal;
    double *oldVal;
    double *tempVal;
    
    feedbackSolutionMDP* optimum;


    //Allocate buffers (to be used to iterate state space and action space).
    int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
    int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

    //Allocate newVal, oldVal, res and optimum action arrays.
    newVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );
    oldVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );
    optimum=new feedbackSolutionMDP();
    optimum->setSize(stateSpace->cardinal());
    // tableau qui pour chaque etat enregistre l'action optimale
    int *action = (int*) calloc( stateSpace->cardinal(), sizeof(int) );
    optimum->setAction(action);


    //Initialization.
    nbIter = 0;
    
    //To iterate the state space, initialize statebuffer by the first state in the state space.
    stateSpace->firstState(statebuffer);
    //Initialize V0.
    for(i = 0 ; i < stateSpace->cardinal() ; i++)
    {
        oldVal[stateSpace->index(statebuffer)] = 0;
        //Move to next state.
        stateSpace->nextState(statebuffer);
    }
    //Start value iterations.
    do
    {
        //increment the number of iterations.
        nbIter++;
        //To iterate the state space, initialize statebuffer by the first state in the state space.
        stateSpace->firstState(statebuffer);
        for(i = 0 ; i < stateSpace->cardinal() ; i++)
        {
            if(type_r == "min")    //If it is a minimization problem, initialize optimal variable with a very big value.
            {
                optimal = MAXDOUBLE;
            }
            else            //If it is a maximization problem, initialize optimal variable with a very small value.
            {
                optimal = MINDOUBLE;
            }
            indiceEtat=stateSpace->index(statebuffer);
            //To iterate the action space, initialize actionbuffer by the first action in the action space.
            actionSpace->firstState(actionbuffer);
            for(a = 0 ; a < actionSpace->cardinal() ; a++)
            {
                indiceAction= actionSpace->index(actionbuffer);
                //Calculate the gain (evaluation function) for the current state at statebuffer.
                //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state);
                //the value of the multiplication evaluateValueState is stored at res.
                //one alternates between oldval and newval to avoid the copy
                res = transitions.at(indiceAction)->evaluateValueState(oldVal, indiceEtat);
                gain = rewards->getEntry(indiceEtat,indiceAction) + res;
                
                //Check if we have a new optimal (min or max depending on the problem type).
                if(((gain < optimal)&&(type_r == "min"))||((gain > optimal)&&(type_r == "max")))
                {
                    //New optimal found, update.
                    optimal = gain;
                    argopt = indiceAction;
                }
                //Move to next action.
                actionSpace->nextState(actionbuffer);
            }
            //Save optimal action and its associated gain value for the current state (statebuffer).
            newVal[indiceEtat] = optimal;
            optimum->setActionIndex(indiceEtat,argopt);
            //Move to next state.
            stateSpace->nextState(statebuffer);
        }
        //Calculate the new distance between Vn+1 and Vn.
        distance = Norm(newVal, oldVal, stateSpace->cardinal());
        // exchange 
        tempVal=oldVal;
        oldVal=newVal;
        newVal=tempVal;
    }while((distance > epsilon) && (nbIter < maxIter));
    printf("# Value Iteration Done %d iterations. Final distance = %e\n", nbIter, distance);


    optimum->setValue(oldVal);
    //Free unnecessary memory blocks.
    free(newVal);
    

    free(actionbuffer);
    free(statebuffer);

    //return the resulting solution.
    return optimum;
}


/* *********************************************************************** */
// A function that works like PI but evaluate the cost of a policy with iteration of the power.
solutionMDP* totalRewardMDP::policyIterationModified(double epsilon, int maxIter, double delta, int maxInIter){
    int i, a;
    int indiceEtat;
    int indiceAction;
    int indiceArgOpt;
    int nbIter =0; //nb iterations totales
    int nbInIter=0; // nb iterations internes

    bool flag=true;
    double distance;
    double gain; // pour calculer la gain d'une action
    double optimal; // garde l'optimal
    double res; // res phase

    double *tmpVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur temporaire
    double *initVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur initiale
    double *newVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );  // fonction de valeur
    double *echVal;

    //Allocate buffers (to be used to iterate state space and action space).
    int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
    int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

    feedbackSolutionMDP *optimum;
    optimum=new feedbackSolutionMDP();
    optimum->setSize(stateSpace->cardinal());
    // tableau qui pour chaque etat enregistre l'action optimale
    int *tmpAction = (int*) calloc( stateSpace->cardinal(), sizeof(int) );
    

    /* boucle generique */
    do {
        //increment the number of iterations.
        nbIter++;
        /*Policy iteration phase. */
        //To iterate the state space, initialize statebuffer by the first state in the state space.
        stateSpace->firstState(statebuffer);
        for(i = 0 ; i < stateSpace->cardinal() ; i++)
        {
            if(type_r == "min")             //If it is a minimization problem, initialize optimal variable with a very big value.
            {
                optimal = MAXDOUBLE;
            }
            else                            //If it is a maximization problem, initialize optimal variable with a very small value.
            {
                optimal = MINDOUBLE;
            }
            indiceEtat=stateSpace->index(statebuffer);
            //To iterate the action space, initialize actionbuffer by the first action in the action space.
            actionSpace->firstState(actionbuffer);
            for(a = 0 ; a < actionSpace->cardinal() ; a++)
            {
                //Calculate the gain (evaluation function) for the current state at statebuffer.
                //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state);
                indiceAction =actionSpace->index(actionbuffer);
                res = transitions.at(indiceAction)->evaluateValueState(initVal,indiceEtat);
                gain = rewards->getEntry(indiceEtat,indiceAction) + res;
                //Check if we have a new optimal (min or max depending on the problem type).
                if(((gain < optimal)&&(type_r == "min"))||((gain > optimal)&&(type_r == "max")))
                {
                    //New optimal found, update.
                    optimal = gain;
                    indiceArgOpt = actionSpace->index(actionbuffer);
                }
                //Move to next action.
                actionSpace->nextState(actionbuffer);
            }
            //Save optimal action and its associated gain value for the current state (statebuffer).
            tmpVal[indiceEtat] = optimal;
            tmpAction[indiceEtat]=indiceArgOpt;
            //Move to next state.
            stateSpace->nextState(statebuffer);
        }
        //Calculate the new distance between Vn+1 and Vn.
        distance = Norm(tmpVal,initVal, stateSpace->cardinal());
        //Optimality check.
        if(distance < epsilon)
        {
            //Optimal found.
            flag = false;
        }
        else
        //Optimal not found, run value iteration to improve.
        {
            nbInIter=0;
            //Value iteration phase.
            do
            {
                nbInIter++;
                //To iterate the state space, initialize statebuffer by the first state in the state space.
                stateSpace->firstState(statebuffer);
                for(i = 0; i < stateSpace->cardinal() ; i++)
                {
                    indiceEtat=stateSpace->index(statebuffer);
                    indiceAction=tmpAction[indiceEtat];
                    res = transitions.at(indiceAction)->evaluateValueState(tmpVal,indiceEtat);
                    newVal[indiceEtat] = rewards->getEntry(indiceEtat,indiceAction) + res;
                    //Move to next state.
                    stateSpace->nextState(statebuffer);
                }
                distance = Norm(newVal, tmpVal, stateSpace->cardinal());
                // echange des pointeurs des vecteurs pour recommencer en 1.
                echVal=newVal;
                newVal=tmpVal;
                tmpVal=echVal;
            }while ( (distance > delta ) && ( nbInIter < maxInIter) );
            // fin du calcul de la valeur
            // echange des espaces pointes
            echVal=initVal;
            initVal=tmpVal;
            tmpVal=echVal;
        } // fin de la phase iteration valeur
    // fin du step
    }while((flag) && (nbIter < maxIter));
    printf("# value Iteration Modified Done %d iterations. Final distance = %e\n", nbIter, distance);

    //
    optimum->setValue(tmpVal);
    optimum->setAction(tmpAction);

    free(newVal);
    free(initVal);
    free(actionbuffer);
    free(statebuffer);

    return optimum;
}


  
