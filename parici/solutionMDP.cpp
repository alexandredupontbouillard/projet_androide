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

#include <list>
#include <vector>
#include <iostream>
#include <sstream>

#include "header/solutionMDP.h"

using namespace std;

//This .cpp file is for implementing functions of the abstract solutionMDP class.

solutionMDP::solutionMDP(){
  _size=0;
}

solutionMDP::~solutionMDP(){
}

void solutionMDP::setSize(int s){
    _size=s;
}

// for printing
void solutionMDP::writeSolution(){
 cout <<"#Print solution of an MDP problem "<< endl;
 cout <<"#Size of the state space : "<<_size<<"\n"<< endl;
}


