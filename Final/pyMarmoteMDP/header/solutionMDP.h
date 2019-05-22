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

Copyright 2018 EMmanuel Hyon, Alain Jean-Marie Abood Mourad*/

#ifndef SOLUTIONMDP_H
#define SOLUTIONMDP_H


/**
 * @brief Class solutionMDP: implementation of an abstract solutionMDP class.
 * 
 * This class is to be inherited by other classes.
 * @author Abood Mourad, and E. Hyon, lip6 
 * @date 18 jan 2018
 * @version 2
 *
 */
class solutionMDP
{
public:

  /**
  * @brief Constructor to create solutionMDP object.
  * @author Hyon
  * @version 1
  * @date feb 2018
  */
  solutionMDP();

  /** 
  * @brief the destructor of the object solutionMDP
  * @author EH
  * @version 1
  * @date feb 2018
  */
  virtual ~solutionMDP();


  /**
  * @brief A function to print the solution object.
  * @author Abood Mourad.
  * @date feb 2018
  * @return none.
  */
   virtual void writeSolution();
   
   
  /**
  * @brief setter of the size.
  * @author EH
  * @date mar 2018
  * @return none.
  */ 
  void setSize(int s);
     

protected:
  int _size;                         /**< _size of the solution */
};

#endif // SOLUTIONMDP_H




