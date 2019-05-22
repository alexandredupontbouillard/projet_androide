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

Copyright 2018 EMmanuel Hyon, Alain Jean-Marie*/


#include <malloc.h>
#include "header/feedbackSolutionMDP.h"

using namespace std;

//This .cpp file is for implementing functions of the feedbackSolutionMDP class.


// constructeur qui initialise deux vecteurs  null
feedbackSolutionMDP::feedbackSolutionMDP(){
    action=NULL;
    value=NULL;
}

// destructor
feedbackSolutionMDP::~feedbackSolutionMDP(){
    if (action != NULL)
        free(action);
    if (value != NULL)
        free(value);
}


//setter action
void feedbackSolutionMDP::setAction(int* a){
   action=a;
}

//setter action
void feedbackSolutionMDP::setValue(double* t){
   value=t;
}

//setter one action (optimal action)
void feedbackSolutionMDP::setActionIndex(int indice,int actionsimple){
  action[indice]=actionsimple;
}

//getter optimal action at index
int feedbackSolutionMDP::getActionIndex(int indice){
  return action[indice];
}

//getter optimal value at index
double feedbackSolutionMDP::getValueIndex(int indice){
  return value[indice];
}

//A function to print the solution object.
void feedbackSolutionMDP::writeSolution()
{
	int i;
	solutionMDP::writeSolution();

	printf("#############################################\n");
	printf("# Solution of the entered problem model:\n");
	printf("# - column 1: index of the state\n");
	printf("# - column 2: Value function \n");
	printf("# - column 3: Optimal action\n");
	printf("#\n");
	for (i = 0 ; i < _size ; i++)
	{
		printf("%2d \t%10.6f %3d", i, value[i], action[i]);
		printf("\n");
	}
	printf("#############################################\n");
}
