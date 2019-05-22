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


/**
   	* @brief library of linear algebra functions
   	* @author Hyon
	* @version 1
	* @date jan 2018
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#ifndef ALGLIN_H
#define ALGLIN_H

extern void  verifVd(double *,int);
extern void  verifVi(int *,int);
extern void  verifMat(double **,int);
extern void  verifMat2D(double **,int,int);
extern void  verifMatP(double **,int);
extern double *produitMatVect(double **,double *,int);
extern double **produitMatMat(double **,double **,int);
extern double prodScalV(double *,int);

/**
  * @brief function to calculate the norm 1 of the difference of two vectors of doubles.
  * @author Hyon
	* @param double* vector 1
	* @param double* vector 2
  * @param size of the vectors
  * @return the difference
*/
extern double Norm( double *, double *, int);

/**
  * @brief function to calculate the span of the difference of two vectors of doubles.
  * @author Hyon
	* @param double* vector 1
	* @param double* vector 2
  * @param size of the vectors
  * @return the difference of the maximum gap minus the minimum gap
*/
extern double Span( double *, double *, int);
/**
  * @brief function to calculate the span of the difference of two vectors of doubles and to get the .
  * @author Hyon
	* @param double* vector 1
	* @param double* vector 2
  * @param size of the vectors
  * @return the difference of the maximum gap minus the minimum gap
*/
extern double SpanRecup( double *U, double *V, int Size, double *max, double *min);

extern void Transpose(double **,int);
extern double *ResolutionSysLin(double **,double *, int);
extern void Inversion(double **,double **, int);
extern double *CalculPI(double **,int);
extern double *CalculPIChaine(double **,int);

#endif // _ALGLIN_H
