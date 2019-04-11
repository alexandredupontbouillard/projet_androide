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

/**
 * @brief Class genericMDP source
 * @author Hyon, lip6
 * @date 18 jan 2018
 * @version 3
 *
 * This .cpp file is for implementing functions of the genericMDP class
 */

#include "header/genericMDP.h"


#include <iostream>
#include <sstream>
#include <malloc.h>

using namespace std;

// ****
// Constructor to create a genericMDP object.
// ce constructeur prend une matrice de reward de taille = taille espace action * taille espace etat cela donne rews(a,x) qui donne le gain r(x,a) dans le modele
genericMDP::genericMDP(string t, string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, sparseMatrix* rews)
{
	type_t = t;
	type_r = r;
	stateSpace = states;
	actionSpace = actions;

	for(int a = 0 ; a < actionSpace->cardinal() ; a++)
	{
      //la methode push_back ajoute l'element a la fin du vecteur
	    transitions.push_back(trans.at(a));
	}
	rewards = rews;
}

// ****
//Second Constructor to create a genericMDP object.
genericMDP::genericMDP(string t, string r,  marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, vector<sparseMatrix*> rews)
{
	type_t = t;
	type_r = r;
	stateSpace = states;
	actionSpace = actions;
	for(int a = 0 ; a < actionSpace->cardinal() ; a++)
	{
	    transitions.push_back(trans.at(a));
	}
	rewards=cost_perStage(rews);
}

// ****
//Destructor for a genericMDP object.
genericMDP::~genericMDP()
{ }

// ****
//Override the virtual function of the abstract class (override the printing function).
void genericMDP::writeMDP()
{
    printf("#############################################\n");
    printf("write MDP\n");
    printf("MDP type (discrete,continuous): %s \n" , type_t.c_str());
    printf("MDP rule (min,max): %s \n",type_r.c_str());
    printf("#############################################\n");
    printf("State space size: %d\n" , (int)stateSpace->cardinal());
    printf("Action space size: %d\n" , (int)actionSpace->cardinal());
    printf("State  dimension: %d\n" , (int)stateSpace->totNbDims());
    printf("Action dimension: %d\n" , (int)actionSpace->totNbDims());
    printf("#############################################\n");
    printf("Transition matrix per action:\n");
    int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
    actionSpace->firstState(actionbuffer);
    for(int k = 0 ; k < actionSpace->cardinal() ; k++)
    {
		printf("action: %d\n" , actionSpace->index(actionbuffer));
		for(int i = 0 ; i < stateSpace->cardinal() ; i++)
		{
			for(int j = 0 ; j < stateSpace->cardinal() ; j++)
			{
		 		printf("%f	", transitions.at(actionSpace->index(actionbuffer))->getEntry(i,j));
			}
			printf("\n");
  		}
		actionSpace->nextState(actionbuffer);
    }
    printf("#############################################\n");
    printf("Reward Matrix (action, state):\n");
    for(int i = 0 ; i < stateSpace->cardinal() ; i++)
    {
		for(int j = 0 ; j < actionSpace->cardinal() ; j++)
		{
			printf("%f	", rewards->getEntry(i,j));
		}
		printf("\n");
    }
    printf("#############################################\n");
    printf("#############################################\n");
    free(actionbuffer);
}

void genericMDP::clearRew(){
  if (NULL != rewards){
     delete rewards;
		 rewards=NULL;
  }
}

// ****
//A function to calculate costs per stage.
sparseMatrix* genericMDP::cost_perStage(vector<sparseMatrix*> rews)
{
	sparseMatrix* rewardsTemporaire = new sparseMatrix(stateSpace->cardinal(),actionSpace->cardinal());
	//Allocate buffers (to be used to iterate state space and action space).
	int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) ); // to scan the action Space
	int* statebuffer1 = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );  // to scan the initial state space
	int* statebuffer2 = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );  // to scan the final state space
	double res = 0.0;
	//To iterate the action space, initialize actionbuffer by the first action in the action space.
	actionSpace->firstState(actionbuffer);
	for(int k = 0 ; k < actionSpace->cardinal() ; k++)
	{
		//To iterate the state space, initialize statebuffer by the first state in the state space.
		stateSpace->firstState(statebuffer1);
		for(int i = 0 ; i < stateSpace->cardinal() ; i++)
		{
			res = 0.0;
			//To iterate the state space, initialize statebuffer by the first state in the state space.
			stateSpace->firstState(statebuffer2);
			for(int j = 0 ; j < stateSpace->cardinal() ; j++)
			{
                // compute the reward associated with the transition
                res += transitions.at(actionSpace->index(actionbuffer))->getEntry(stateSpace->index(statebuffer1),stateSpace->index(statebuffer2)) * rews.at(actionSpace->index(actionbuffer))->getEntry(stateSpace->index(statebuffer1),stateSpace->index(statebuffer2));
				//Move to next state.
				stateSpace->nextState(statebuffer2);
			}
			//Set the accumulated cost per stage to the rewards matrix.
			rewardsTemporaire->setEntry(stateSpace->index(statebuffer1),actionSpace->index(actionbuffer),res);
			//Move to next state.
			stateSpace->nextState(statebuffer1);
		}
		//Move to next action.
		actionSpace->nextState(actionbuffer);
	}
    //Free unnecessary memory blocks.
    free(actionbuffer);
    free(statebuffer1);
    free(statebuffer2);

    return rewardsTemporaire;
}
