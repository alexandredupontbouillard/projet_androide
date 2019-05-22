/* Marmote is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Marmote is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Marmote. If not, see <http://www.gnu.org/licenses/>.

Copyright 2015 Alain Jean-Marie, Jean-Michel Fourneau, Jean-Marc Vincent, Issam Rabhi */

#include "transitionStructure.h"
#include <sstream>
#include <string>
#include <iomanip> 
#include <stdlib.h>
//#include <RcppArmadillo.h>      	// for Armadillo as well as Rcpp 
//#include <RInside.h>                    // for the embedded R via RInside

// using namespace Rcpp;

discreteDistribution* transitionStructure::evaluateMeasure(discreteDistribution* d) {

  double *measure;
  discreteDistribution *res;

  measure = (double*)malloc( _origSize * sizeof(double) );
  evaluateMeasure( d->probas(), measure );
  res = new discreteDistribution( _origSize, d->values(), measure );
  free( measure );
  
  return res;
}

bool transitionStructure::readEntry(FILE *input)
{
  double	val;
  int		row;
  int		col;

  if ( 3 == fscanf( input, "%d %d %lf", &row, &col, &val ) ) {
    return setEntry( row, col, val );
  }
  else {
    return false;
  }
}

// Default implementation of the consolidation procedure. Assumes that the
// matrix is full as a worst-case.
// Is destined to be reimplemented for specific structures.
int transitionStructure::consolidate(int i, int *destinations, double *values)
{
  int res = _origSize;

  for ( int j = 0; j < _destSize; j++ ) {
    destinations[j] = j;
    values[j] = getEntry(i,j);
  }

  return res;
}

// Default implementation of the write procedure. Assumes that the
// matrix is full as a worst case.
// Is destined to be reimplemented for specific structures.
void transitionStructure::write(FILE *out, std::string format)
{

  if ( NULL == out ) {
    return;
  }

  // assert: out is not NULL

  if (format == "XBORNE"){
    fprintf( out, "%s", toString("R").c_str() );
    // need to create modelName.sz (and modelName.cd)
  }
  else if ( ( format == "MARCA" ) || ( format == "MatrixMarket-sparse" ) ) {
    // both MARCA and Matrix Market are "1-based"

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii+1, jj+1, pro ); }

    // enumerate all states
    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      for ( int d = 0; d < _destSize; d++ ) {
        FORMAT( i, d, getEntry(i,d) );
      }
    }
  }
  else if ( format == "Ers" ) {

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "%10d %10d %12e\n", ii, jj, pro ); }

    for ( int i=0; i < _origSize; i++ ) {
      for ( int j=0; j < _destSize; j++ ) {
        FORMAT( i, j, getEntry(i,j) );
      }
    }
  }
  else if ( format == "Maple" ) {
    // Maple matrices are "1-based".

#undef FORMAT
#define FORMAT( ii, jj, pro ) { fprintf( out, "(%d,%d)=%12e,\n", ii+1, jj+1, pro ); }
#define FORMAT2( ii, jj, pro ) { fprintf( out, "(%d,%d)=%12e\n", ii+1, jj+1, pro ); }

// error#2: initialize nbTrans;

int nbTrans=10;

    for ( int i = 0; i < _origSize; i++ ) {
      // handle state
      if ( i < _origSize-1 ) {
        for ( int d = 0; d < _destSize; d++ ) {
          FORMAT( i, d, getEntry(i,d) );
        }
      }
      else {
        for ( int d = 0; d < _destSize-1; d++ ) {
          FORMAT( i, d, getEntry(i,d) );
        }
        // avoid comma at the end of last record
        FORMAT2( i, (int)_destSize-1, getEntry(i,nbTrans-1) );
      }
    }
  }
  else if ( format == "MatrixMarket-full" ) {
    // format is "column-major" (!)
    for ( int j=0; j<_destSize; j++ ) {
      for ( int i=0; i<_origSize; i++ ) {
        double val = getEntry( i, j );
        fprintf( out, " %12e", val );
      }
      fprintf( out, "\n");
    }
  }
  else if ( format == "Full" ) {
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_destSize; j++ ) {
        double val = getEntry( i, j );
        fprintf( out, " %10.8f", val );
      }
      fprintf( out, "\n");
    }
  }
  else if ( format == "R" ) {
    fprintf( out, "%s", toString("R").c_str() );
  }
 else if ( format == "SCILAB" ) {
   fprintf( out, "%s", toString("SCILAB").c_str() );

    // Scilab full format. Not complete.
 /*   fprintf( out, "full( [ ");
    double val;
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_destSize-1; j++ ) {
        val = getEntry( i, j );
        fprintf( out, " %10.8f,", val );
      }
      if ( i < _origSize-1 ) {
        val = getEntry( i, _origSize-1 );
        fprintf( out, " %10.8f;", val );
      }
      else {
        fprintf( out, " %10.8f]", val );
      }
    }
*/
  }
  else if ( format == "Matlab" ) {
    // Matlab full format
    // Not complete
    // Finally print array of values
    fprintf( out, "theValues = [");
    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_destSize; j++ ) {
        fprintf( out, " %12f", getEntry(i,j) );
      }
    }
    fprintf( out, "]';\n" );

  }
  else {
    fprintf( stderr, "Warning in transitionStructure::write(): format %s not recognized. Ignored.",
             format.c_str() );
  }

}

std::string transitionStructure::toString( std::string format )
{
  std::string res;
  
  if ( format == "R" ) {
    
    std::string mat="matrix(c("; // string to store the matrix

    for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_destSize; j++ ) {
        double val = getEntry(i,j);
        std::ostringstream tmp;
        tmp.precision(10);
        tmp << val;
        // end of lines
        if ( (j==_destSize-1)&&(i==_origSize-1) ) {
          mat.append(tmp.str());
        }
        else {
          mat.append(tmp.str()+",");
        }
      }
    }
    std::stringstream convert;
    convert << _destSize;
    res = mat + "), nrow=" + convert.str() + ", byrow=TRUE)";
  }
  else if ( format == "SCILAB" ) {

	std::string mat="mcA=[";  // string to store the matrix

     for ( int i=0; i<_origSize; i++ ) {
      for ( int j=0; j<_destSize; j++ ) {
        double val = getEntry(i,j);
        std::ostringstream tmp;
        tmp.precision(10);
        tmp << val;
        // end of lines
        if ( (j==_destSize-1)&&(i==_origSize-1) ) {
          mat.append(tmp.str());
        }
        else {
          mat.append(tmp.str()+" ");
        }
      }
      if(i!=_origSize-1) mat.append(";");
     }
     std::stringstream convert;
     convert << _destSize;
     res = mat + "]";

  }
  else if ( format == "XBORNE" ) {

    std::string mat = "";

    for ( int i=0; i<_origSize; i++ ) {
      if ( getNbElts(i) > 0 ) {
        std::ostringstream line;
        line << i;
        for ( int j=0; j<_destSize; j++ ) {
          double val = getEntry(i,j);
          if ( val != 0.0 ) {
            line << " ";
            line.width(12);
            line << j;
            line << " ";
            line.width(12);
            line.precision(10);
            line << val;
          }
        }
        mat.append( line.str() +"\n");
      }
    }
    res = mat;
  }
  else {
    fprintf( stderr, "Warning in transitionStructure::toString(): format %s not recognized. Ignored.",
             format.c_str() );
    res = "";
  }

  return res;
}
