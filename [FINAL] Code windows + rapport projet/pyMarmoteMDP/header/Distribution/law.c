/************************************************************************/
/* Module containing functions for manipulating distribution laws.	*/
/* Provides:								*/
/*	void	Write_Law( FILE *Stream, Law_Desc Loi, int Mode );	*/
/*	double	Mean( Law_Desc Loi );					*/
/*	double	Rate( Law_Desc Loi );					*/
/*	double	Moment( Law_Desc Loi, int order );			*/
/*	double	Laplace( Law_Desc Loi, double s );			*/
/*	double	DLaplace( Law_Desc Loi, double s );			*/
/*	Law_Desc New_Law( void );					*/
/*	Law_Desc Copy_Law( law_Desc Loi );				*/
/*	Law_Desc Rescale_Law( Law_Desc Loi, double Factor );		*/
/*      Law_Desc Parse_Law( char Nom, Liste_Reel Params, 		*/
/*			Liste_Reel Params2, Law_List The_SubLaws );	*/
/*	Law_Desc Parse_Law_From_Args ( char**, int, int, Law_Desc );	*/
/*	Law_List Append_Law( Law_List, Law_Desc New_Law );		*/
/*	void Free_Law( Law_Desc The_Law );				*/
/************************************************************************/
#define	TOOLNAME	"LawPackage"

#include <stdlib.h>
#include <math.h>
#include "law.h"
#include "messages.h"
#define M_SQRT1_2   0.70710678118654752440  /* 1/sqrt(2) */

#define SCAN(x,y)	for(x=y;x!=NULL;x=x->Next)
#define EXPAND(ptr,nb)	ptr = (double*) realloc( ptr, (size_t)(nb)*sizeof(double) )

/* external declarations */
double erfc(double x);
double rint(double x);
double lgamma(double x);

/* Constructor. The parameter table is initialized with room for one. */
Law_Desc New_Law(void) {

  Law_Desc The_Law;

  The_Law = (Law_Desc) malloc( sizeof( struct LAW_DESC ) );
  The_Law->Parameters = (double*) malloc( sizeof(double) );
  The_Law->Law_Parameters = (Law_Desc*) NULL;

  return The_Law;
}

/* Procedure to write the law. Accepts various styles.			*/
void Write_Law( FILE *out, 
		Law_Desc The_Law, 
		int Mode )
{
  int	Nb_Params;
  int	i;

  switch( Mode ) {
  case PNED_PRINT_MODE:
    fprintf( out, "{ name " );
    switch( The_Law->Name ) {
    case 'A':
      /* Convention: shape = alpha, rate = beta */
      fprintf( out, "Gamma shape %8.4f rate %8.4f ", 
	       The_Law->Parameters[0], The_Law->Parameters[1] );
      break;
    case 'G':
      fprintf( out, "Geometric proba %8.4f ", The_Law->Parameters[0] );
      break;
    case 'I': 
      fprintf( out, "constant infinite value " );
      break;
    case 'W':
      fprintf( out, "weibull a %8.4f b %8.4f", 
	       The_Law->Parameters[0], The_Law->Parameters[1] );
      break;
    case 'N':
      fprintf( out, "normal mean %8.4f stddev %8.4f", 
	       The_Law->Parameters[0], The_Law->Parameters[1] );
      break;
    case 'L':
      fprintf( out, "lognormal mean %8.4f stddev %8.4f", 
	       The_Law->Parameters[0], The_Law->Parameters[1] );
      break;
    case 'p':
      fprintf( out, "pareto alpha %8.4f beta %8.4f ", 
	       The_Law->Parameters[0], The_Law->Parameters[1] );
      break;
    case 'R':
      fprintf( out, "erlang phases %4.0f mean %8.4f ", 
	       The_Law->Parameters[0], The_Law->Parameters[1] );
    break;
    case 'O':
      fprintf( out, "%c %8.4f ", The_Law->Name, The_Law->Parameters[0] );
      fprintf( out, "%8.4f %8.4f ",
	       The_Law->Parameters[1],
	       The_Law->Parameters[2] );
      break;
    case 'H':
      fprintf( out, "hyperexpo means { ");
      for ( i=1; i<The_Law->Parameters[0]; i+=2 ) {
	fprintf( out, "%8.4f ", The_Law->Parameters[i] );
      }
      fprintf( out, "} probas { ");
      for ( i=2; i<The_Law->Parameters[0]; i+=2 ) {
	fprintf( out, "%8.4f ", The_Law->Parameters[i] );
      }
      fprintf( out, "} ");
      break;
    case 'D':
      fprintf( out, "discrete values { ");
      Nb_Params = (int)rint(The_Law->Parameters[0]);
      for ( i=1; i<= Nb_Params; i++ )
	fprintf( out, "%8.4f ", The_Law->Parameters[i] );
      fprintf( out, "} probas { ");
      for ( i=1; i<= Nb_Params; i++ )
	fprintf( out, "%8.4f ", The_Law->Parameters[i+Nb_Params] );
      fprintf( out, "} ");
      break;
    case 'P':
      fprintf( out, "IPP [ " );
      fprintf( out, "%f, ", The_Law->Parameters[0] );
      Write_Law( out, The_Law->Law_Parameters[0], Mode );
      fprintf( out, ", ");
      Write_Law( out, The_Law->Law_Parameters[1], Mode );
      fprintf( out, " ] ");
      break;
    default:
      fprintf( stderr, "%s: error (Print): unexpected law id (%c)\n",
	       TOOLNAME, The_Law->Name );
    }
    fprintf( out, "}" );
    break;

  default:
    switch( The_Law->Name ) {
    case 'I':
      fprintf( out, "%c ", The_Law->Name );
      break;
    case 'g':
      fprintf( out, "%c %8.4f ", The_Law->Name, The_Law->Parameters[0] );
      break;
    case 'A':
    case 'W':
    case 'N':
    case 'L':
    case 'p':
      fprintf( out, "%c %8.4f %8.4f ", 
	       The_Law->Name,
	       The_Law->Parameters[0],
	       The_Law->Parameters[1] );
      break;
    case 'R':
      fprintf( out, "%c %4.0f %8.4f ", The_Law->Name,
	       The_Law->Parameters[0],
	       The_Law->Parameters[1] );
      break;
    case 'O':
      fprintf( out, "%c %8.4f ", The_Law->Name, The_Law->Parameters[0] );
      fprintf( out, "%8.4f %8.4f ",
	       The_Law->Parameters[1],
	       The_Law->Parameters[2] );
      break;
    case 'H':
      fprintf( out, "%c ", The_Law->Name );
      for ( i=0; i<The_Law->Parameters[0]; i+=2 ) {
	fprintf( out, "%8.4f %8.4f",
		 The_Law->Parameters[i+1],
		 The_Law->Parameters[i+2] );
      }
      break;
    case 'D':
      fprintf( out, "%c [ ", The_Law->Name );
      Nb_Params = (int)rint(The_Law->Parameters[0]);
      for ( i=1; i<= Nb_Params; i++ )
	fprintf( out, "%8.4f ", The_Law->Parameters[i] );
      fprintf( out, "] [ ");
      for ( i=1; i<= Nb_Params; i++ )
	fprintf( out, "%8.4f ", The_Law->Parameters[i+Nb_Params] );
      fprintf( out, "] ");
      break;
    case 'P':
      fprintf( out, "IPP [ " );
      fprintf( out, "%f, ", The_Law->Parameters[0] );
      Write_Law( out, The_Law->Law_Parameters[0], Mode );
      fprintf( out, ", ");
      Write_Law( out, The_Law->Law_Parameters[1], Mode );
      fprintf( out, " ] ");
      break;
    default:
      fprintf( stderr, "%s: error (Print): unexpected law id (%c)\n",
	       TOOLNAME, The_Law->Name );
    }
  }
}

double Mean(  Law_Desc The_Law )
{	
  double	The_Mean;
  int		i;
  int		Nb_Params;
  double	Alpha;
  double       	Mean_On;
  double       	Mean_Off;

  switch( The_Law->Name ) {
  case 'I':
    The_Mean = INFINITE_DURATION;
    break;
  case 'O':
  case 'N':
  case 'L':
    The_Mean = The_Law->Parameters[0];
    break;
  case 'A':
    The_Mean = The_Law->Parameters[0] / The_Law->Parameters[1];
    break;
  case 'l':
    The_Mean = 2 * The_Law->Parameters[0];
    break;
  case 'G':
    The_Mean = 1.0 / ( 1.0 - The_Law->Parameters[0] );
    break;
  case 'R':
    The_Mean = The_Law->Parameters[1];
    break;
  case 'p':
    if ( The_Law->Parameters[1] > 1.0 )
      The_Mean = The_Law->Parameters[0]/(The_Law->Parameters[1] - 1.0);
    else
      The_Mean = 1.0/0.0;
    break;
  case 'W':
    Alpha = 1.0 / The_Law->Parameters[1];
    The_Mean = Alpha * exp( lgamma( Alpha ) 
			    - Alpha * log( The_Law->Parameters[0] ) );
    break;
  case 'H':
    The_Mean = 0.0;
    for ( i=0; (i<The_Law->Parameters[0]); i+=2 )
      The_Mean += The_Law->Parameters[i+1]*The_Law->Parameters[i+2];
    break;
  case 'P':
    Mean_On = Mean( The_Law->Law_Parameters[0] );
    Mean_Off = Mean( The_Law->Law_Parameters[1] );
    The_Mean = The_Law->Parameters[0] * Mean_On / ( Mean_On + Mean_Off );
    break;
  default:
    fprintf( stderr, "%s: warning: cannot compute mean of law `%c`.",
	     TOOLNAME, The_Law->Name );
    fprintf( stderr, " Zero assumed.\n");
    The_Mean = 0.0;
  }

  return	The_Mean;
}

double
Rate( Law_Desc The_Law )
{
  double	The_Rate;
  double	The_Mean;

  switch( The_Law->Name ) {
  case 'I':
    The_Rate = 0.0;
    break;
  default:
    The_Mean = Mean( The_Law );
    if ( The_Mean == 0.0 )
      The_Rate = INFINITE_RATE;
    else
      The_Rate = 1.0 / The_Mean;
  }

  return The_Rate;
}

double 
Moment(       Law_Desc		The_Law,
	      int		order )
{	
  double	The_Moment = 0.0;
  double	s2;
  double	m;
  int		i;
  int		Nb_Params;
  double	*Param;

  Param = The_Law->Parameters;

  if ( order == 0 )
    return 1.0;
  else if ( order == 1 )
    return Mean( The_Law );

  /* from here on, order >= 2 */
  switch( The_Law->Name ) {
  case 'A':
    The_Moment = Param[0] / Param[1];
    for ( i = 1; i < order; i++ ) {
      The_Moment = The_Moment * ( Param[0] + i ) / Param[1];
    }
    break;
  case 'G':
    if ( order == 2 )
      The_Moment = ( 1.0 + Param[0] ) / ( 1.0 - Param[0] );
    else {
      fprintf( stderr, "%s: warning: cannot compute moment of order > 1",
	       TOOLNAME );
      fprintf( stderr, "for law `%c`. 0.0 assumed.\n", The_Law->Name );
    }
    break;
  case 'I':
    The_Moment = INFINITE_DURATION;
    break;
  case 'R':
    The_Moment = pow( Param[1], order );
    for( i=1; i<=order-1; i++ )
      The_Moment *= (Param[0]+(double)i)/Param[0];
    break;
  case 'H':
    The_Moment = 0.0;
    for ( i=0; (i<Param[0]); i+=2 )
      The_Moment += Param[i+1] * pow( Param[i+2], order );
    for( i=1; i<=order; i++ )
      The_Moment *= (double)i;
    break;
  case 'N':
    switch( order ) {
    case 2:
      The_Moment = Param[0]*Param[0] + Param[1]*Param[1];
      break;
    case 3:
      The_Moment = Param[0] *( Param[0]*Param[0] + 3.0*Param[1]*Param[1] );
      break;
    case 4:
      The_Moment = pow( Param[0], 4.0 ) 
	+ 6.0 * pow( Param[0]*Param[1], 2.0 )
	+ 3.0 * pow( Param[1], 4.0 );
      break;
    default:
      fprintf( stderr, "%s: warning: cannot compute %d-th moment of law `%c`.",
	       TOOLNAME, order, The_Law->Name );
      fprintf( stderr, " Zero assumed.\n");
    }
    break;
  case 'L':
    s2 = log( 1.0 + ( Param[1] / Param[0] )*( Param[1] / Param[0] ) );
    m = log( Param[0] ) - s2 / 2;
    The_Moment = exp( order*m + order*order*s2/2.0 );
    break;
  default:
    fprintf( stderr, "%s: warning: cannot compute %d-th moment of law `%c`.",
	     TOOLNAME, order, The_Law->Name );
    fprintf( stderr, " Zero assumed.\n");
    The_Moment = 0.0;
  }

  return	The_Moment;
}

boolean hasMoment(
		  Law_Desc	The_Law,
		  int		order )
{	
  double	*Param;

  Param = The_Law->Parameters;

  if ( order == 0 )
    return TRUE;

  /* from here on, order >= 1 */
  switch( The_Law->Name ) {
  case 'A':
  case 'R':
  case 'H':
  case 'N':
    return TRUE;
  case 'G':
    return ( Param[0] < 1.0 );
  case 'p':
    return ( Param[1] > order );
  case 'I': 
    return FALSE;
    break;
  default:
    fprintf( stderr, "%s: warning: cannot determine whether law `%c` has %d-th moment.",
	     TOOLNAME, The_Law->Name, order );
    fprintf( stderr, " True assumed.\n");
    return TRUE;
  }
}

double Laplace( 
	       Law_Desc	The_Law,
	       double	s )
{	
  double	The_Value;
  int		i;
  int		Nb_Params;
  double	*Param;

  Param = The_Law->Parameters;

  switch( The_Law->Name ) {
  case 'A':
    The_Value = pow( Param[1] / ( Param[1] + s ), Param[0] );
    break;
  case 'G':
    The_Value = ( 1.0 - Param[0] ) * exp( -s )
      / ( 1.0 - Param[0] * exp( -s ) );
    break;
  case 'I':
    The_Value = ( s == 0.0 ? 1.0 : 0.0 );
    break;
  case 'R':
    The_Value = Param[0] / ( s*Param[1] + Param[0] );
    The_Value = exp( Param[0] * log( The_Value ) );
    break;
  case 'H':
    The_Value =  0.0;
    for( i=0; i<Param[0]; i+=2 )
      The_Value += Param[i+1] / ( 1.0 + s * Param[i+2] );
    break;
  default:
    fprintf( stderr, "%s: warning: cannot compute Laplace transform ",
	     TOOLNAME );
    fprintf( stderr, "of law `%c`. 1.0 assumed.\n", 
	     The_Law->Name );
    The_Value = 1.0;
  }

  return	The_Value;
}

double DLaplace( 
	       Law_Desc	The_Law,
	       double	s )
{	
  double	The_Value;
  int		i;
  int		Nb_Params;
  double	*Param;

  Param = The_Law->Parameters;

  switch( The_Law->Name ) {
  case 'A':
    The_Value = Param[0] / Param[1] 
      * pow( Param[1] / ( Param[1] + s ), Param[0] + 1 );
    break;
  case 'G':
    The_Value = - ( 1.0 - Param[0] ) * exp( -s )
      / ( 1.0 - Param[0] * exp( -s ) )
      / ( 1.0 - Param[0] * exp( -s ) );
    break;
  case 'R':
    The_Value =  Param[0] / ( s*Param[1] + Param[0] );
    The_Value = - Param[0]*Param[1] * exp( Param[0] * log( The_Value ) )
      / ( s*Param[1] + Param[0] );
    break;
  case 'H':
    The_Value =  0.0;
    for( i=0; i<Param[0]; i+=2 )
      The_Value -= Param[i+1]*Param[i+2]
	/ (1.0 + s * Param[i+2] ) / (1.0 + s * Param[i+2] );
    break;
  default:
    fprintf( stderr, "%s: warning: cannot compute derivative of Laplace transform ",
	     TOOLNAME );
    fprintf( stderr, "of law `%c`. 1.0 assumed.\n", 
	     The_Law->Name );
    The_Value = 1.0;
  }

  return	The_Value;
}

/* The CDF of a v.a. X is the function \bar F(x) = P( X > x ) */
double CDF( Law_Desc	The_Law,
	    double	x )
{	
  double	The_Value;
  int		i;
  int		Nb_Params;
  double	m;
  double	s2;
  double	*Param;

  Param = The_Law->Parameters;

  switch( The_Law->Name ) {
    /* lazy... will do it later ...
       case 'B':
       The_Value = 1.0 - ( 1.0 - exp( -s ) ) * Param[0];
       break;
       case 'G':
       The_Value = ( 1.0 - Param[0] ) * exp( -s )
       / ( 1.0 - Param[0] * exp( -s ) );
       break;
    */
  case 'I':
    The_Value = 1.0;
    break;
    /* laaaaazy
       case 'R':
       The_Value = Param[0] / ( s*Param[1] + Param[0] );
       The_Value = exp( Param[0] * log( The_Value ) );
       break;
    */
  case 'H':
    if ( x < 0.0 )
      The_Value = 1.0;
    else {
      The_Value =  0.0;
      for( i=0; i<Param[0]; i+=2 ) {
	if ( Param[i+2] > 0.0 ) 
	  The_Value += Param[i+1] * exp( - x * Param[i+2] );
      }
    }
    break;
  case 'W':
    if ( x < 0.0 )
      The_Value = 1.0;
    else 
      The_Value = exp( - Param[0] * pow( x, Param[1] ) );
    break;
  case 'p':
    if ( x < 0.0 )
      The_Value = 1.0;
    else
      The_Value = pow( Param[0]/(Param[0]+x), Param[1] );
    break;
  case 'N':
    The_Value = erfc( ( x - Param[0] ) / Param[1] * M_SQRT1_2 ) / 2.0;
    break;
  case 'L':
    if ( x <= 0.0 )
      The_Value = 1.0;
    else {
      s2 = log( 1.0 + ( Param[1] / Param[0] )*( Param[1] / Param[0] ) );
      m = log( Param[0] ) - s2 / 2;
      The_Value = erfc( ( log(x) - m ) / sqrt( 2 * s2 ) ) / 2.0;
    }
    break;
  default:
    fprintf( stderr, "%s: warning: cannot compute CDF ",
	     TOOLNAME );
    fprintf( stderr, "of law `%c`. 0.0 assumed.\n", 
	     The_Law->Name );
    The_Value = 1.0;
  }

  return	The_Value;
}

/* A procedure to rescale, a law. Creates a new one.		*/ 
Law_Desc Rescale_Law( Law_Desc The_Law, double Factor )
{	
  int		i;
  int		Nb_Params;
  Law_Desc	Res;

  Res = New_Law();

  if ( NULL != Res ) {
    Res->Name = The_Law->Name;

    /* a first pass to allocate space */
    switch( The_Law->Name ) {
    case 'G':
      /* No need to alloc space for 1 parameter. */
      break;
    case 'R':
    case 'W':
    case 'N':
    case 'L':
    case 'p':
      EXPAND( Res->Parameters, 2 );
      break;
    case 'O':
      EXPAND( Res->Parameters, 3 );
      break;
    case 'H':
      Nb_Params = (int)rint(The_Law->Parameters[0]);
      EXPAND( Res->Parameters, 1+Nb_Params );
      break;
    case 'P':
      EXPAND( Res->Parameters, 2 );
      Res->Law_Parameters = (Law_Desc*) calloc( 2, sizeof( Law_Desc ) );
      break;
    default:
      fprintf( stderr, "%s: error: unknown law `%c`. No rescaling.\n",
	       TOOLNAME, The_Law->Name );
      free( Res );
      Res = (Law_Desc)NULL;
    }
  }
  else {
    fprintf( stderr, "%s: error: could not allocate memory for law.\n",
	     TOOLNAME );
    Res = (Law_Desc)NULL;
  }

  if ( Res != NULL ) {
    /* second pass to fill in parameters */
    switch( The_Law->Name ) {
    case 'A':
      Res->Parameters[0] = The_Law->Parameters[0];
      Res->Parameters[1] = The_Law->Parameters[1] / Factor;
      break;
    case 'O':
      Res->Parameters[0] = Factor * The_Law->Parameters[0];
      Res->Parameters[1] = Factor * The_Law->Parameters[1];
      Res->Parameters[2] = Factor * The_Law->Parameters[2];
      break;
    case 'g':
    case 'G':
      fprintf( stderr, "%s: warning: geometric distributions cannot",
	       TOOLNAME );
      fprintf( stderr, "be rescaled. No change.\n");
      Res->Parameters[0] = The_Law->Parameters[0];
      break;
    case 'p':
      Res->Parameters[0] = Factor * The_Law->Parameters[0];
      Res->Parameters[1] = The_Law->Parameters[1];
      break;
    case 'R':
      Res->Parameters[0] = The_Law->Parameters[0];
      Res->Parameters[1] = Factor * The_Law->Parameters[1];
      break;
    case 'W':
      Res->Parameters[0] = The_Law->Parameters[0] 
	* pow( Factor, The_Law->Parameters[1] );
      Res->Parameters[1] = The_Law->Parameters[1];
      break;
    case 'H':
      Res->Parameters[0] = The_Law->Parameters[0];
      for ( i=0; (i<The_Law->Parameters[0]); i+=2 ) {
	Res->Parameters[i+1] = The_Law->Parameters[i+1];
	Res->Parameters[i+2] = Factor * The_Law->Parameters[i+2];
      }
      break;
    case 'P':
      Res->Parameters[0] = Factor * The_Law->Parameters[0];
      Res->Law_Parameters[0] = Rescale_Law( The_Law->Law_Parameters[0],
					    Factor );
      Res->Law_Parameters[1] = Rescale_Law( The_Law->Law_Parameters[1],
					    Factor );
      break;
    default:
      fprintf( stderr, "%s: error: don't know how to rescale law `%c`.\n",
	       TOOLNAME, The_Law->Name );
      free( Res );
      Res = (Law_Desc)NULL;
    }
  }

  return Res;
}

/* A procedure to copy a law. Uses "rescale" with factor 1.0!		*/ 
Law_Desc Copy_Law( Law_Desc The_Law )
{	
  Law_Desc res;

  /* Default value, if anything goes wrong. */
  res = (Law_Desc)NULL;

  /* Checking that the law is "registered". This is only to preempt	*/
  /* the warning message of Rescale_Law, which would sound cryptic...	*/
  switch( The_Law->Name ) {
  case 'A':
  case 'G':
  case 'g':
  case 'p':
  case 'R':
  case 'W':
  case 'O':
  case 'H':
  case 'P':
  case 'N':
  case 'L':
    res = Rescale_Law( The_Law, 1.0 );
    break;
  default:
    fprintf( stderr, "%s: error: unknown law `%c`. No copy.\n",
	     TOOLNAME, The_Law->Name );
    free( res );
    res = (Law_Desc)NULL;
  }

  return res;
}


Law_Desc Parse_Law( 
		   char		Nom, 
		   Liste_Reel 	The_Params,
		   Liste_Reel 	The_Params2,
		   Law_List	The_SubLaws )
{	
  int		Nb_Params;
  int		i;
  Law_Desc	The_Law;
  Liste_Reel	R_Scan;
  Law_List	LL_Scan;

  The_Law = New_Law();

  if ( The_Law != NULL )
    switch( The_Law->Name = Nom ) {
    case 'I':
      /* no parameters at all */
      break;
    case 'g':
    case 'G':
      /* Space for one parameter is already available from New_Law() */
      if ( The_Params != NULL ) {
	The_Law->Parameters[0] = The_Params->Val;
	if ( The_Params->Next != NULL )
	  W_EXT( The_Law->Name );
      }
      else {
	W_MIS( The_Law->Name );
	The_Law->Parameters[0] = 0.0;
      }
      break;
    case 'A':
    case 'W':
    case 'L':
    case 'p':
    case 'h':
      /* 'h' is the "old" QNAP way of defining hyperexponentials by	*/
      /* mean and variance or coefficient of variation. Inactive at the	*/
      /* moment.							*/
      EXPAND( The_Law->Parameters, 2 );
      if ( ( The_Params != NULL ) && ( The_Params->Next != NULL ) ) {
	The_Law->Parameters[1] = The_Params->Val;
	The_Law->Parameters[0] = The_Params->Next->Val;
	if ( The_Params->Next->Next != NULL )
	  W_EXT( The_Law->Name );
      }
      else {
	W_MIS( The_Law->Name );
	The_Law->Parameters[1] = The_Law->Parameters[0] = 0.0;
      }
      break;
    case 'H':
      /* first count the parameters */
      Nb_Params = 0;
      SCAN( R_Scan, The_Params )
	Nb_Params++;
      
      if ( Nb_Params % 2 != 0 ) {
	printf( "Warning: odd number of parameters for HEXP. Last ignored.\n");
	Nb_Params--;
      }
      if ( Nb_Params == 0 ) {
	printf( "Error: not enough parameters for HEXP.\n");
      }
      else {
	/* one has to add a cell to store nb of parameters */
	EXPAND( The_Law->Parameters, 1+Nb_Params );
	The_Law->Parameters[0] = Nb_Params;
	/* store the parameters, reversing the list! */
	SCAN( R_Scan, The_Params )
	  The_Law->Parameters[Nb_Params--] = R_Scan->Val;
      }
      break;
    case 'O':
      EXPAND( The_Law->Parameters, 3 );
      if ( ( The_Params != NULL ) 
	   && ( The_Params->Next != NULL )
	   && ( The_Params->Next->Next != NULL ) ) {
	The_Law->Parameters[2] = The_Params->Val;		/* 1/beta */
	The_Law->Parameters[1] = The_Params->Next->Val;		/* 1/alpha */
	The_Law->Parameters[0] = The_Params->Next->Next->Val;	/* 1/lambda */
	if ( The_Params->Next->Next->Next != NULL )
	  W_EXT( The_Law->Name );
      }
      else {
	W_MIS( The_Law->Name );
      }
      break;
    case 'R': /* differe de U dans l'ordre d'empilement des parametres */
      EXPAND( The_Law->Parameters, 2 );
      if ( ( The_Params != NULL ) && ( The_Params->Next != NULL ) ) {
	The_Law->Parameters[0] = The_Params->Val;
	The_Law->Parameters[1] = The_Params->Next->Val;
	if ( The_Params->Next->Next != NULL )
	  W_EXT( The_Law->Name );
      }
      else {
	W_MIS( The_Law->Name );
	The_Law->Parameters[1] = The_Law->Parameters[0] = 0.0;
      }
      break;

    case 'D':
      /* find the length of the (first) parameter list */
      Nb_Params = 0;
      SCAN( R_Scan, The_Params )
	Nb_Params++;
      EXPAND( The_Law->Parameters, 1 + 2*Nb_Params );
      The_Law->Parameters[0] = Nb_Params;

      /* Store values of the discrete distribution */
      i = 0;
      SCAN( R_Scan, The_Params ) {
	i++;
	The_Law->Parameters[i] = R_Scan->Val;
      }
      /* Store probas of the distribution */
      for ( R_Scan = The_Params2; ( i < 2*Nb_Params ) && ( R_Scan != NULL );
	    R_Scan = R_Scan->Next ) {
	i++;
	The_Law->Parameters[i] = R_Scan->Val;
      }
      if ( R_Scan != NULL ) {
	printf("Warning! Probability list too long (ignored)\n" );
      }
      else if ( i < 2*Nb_Params ) {
	printf("Warning! Probability list too short. 0.0 assumed\n" );
	for( ; i <= 2*Nb_Params; i++ )
	  The_Law->Parameters[i] = 0.0;
      }

    break;
    case 'P':
      The_Law->Parameters[0] = The_Params->Val;
      Nb_Params = 0;
      SCAN( LL_Scan, The_SubLaws )
	Nb_Params++;
      The_Law->Law_Parameters = (Law_Desc*) calloc( (unsigned)Nb_Params,
						   sizeof( Law_Desc ) );
      i = 0;
      SCAN( LL_Scan, The_SubLaws ) {
	The_Law->Law_Parameters[i] = LL_Scan->Val;
	i++;
      }
      break;

    default:
      printf("Warning! Law '%c' ID not recognized\n", The_Law->Name);
    }

  Free_Real_List( The_Params );
  Free_Real_List( The_Params2 );

  return The_Law;
}

/************************************************************************/
/* Procedure to parse a law from the arguments passed to an application	*/
/* on the command line. The arguments are 
   - argv	the standard array of strings
   - index 	the index in this array at which to start,
   - argc	the standard total count of arguments
   - The_Law	a pointer to an _already allocated_ law description
   The return value is the index just after the last token accepted as
   part of the description of the law. Hence, the typical usage will be,
   i being the current index in the argv:
   i = Parse_Law_From_Args( argv, i, argc, The_Trans_Law );		*/
/************************************************************************/
int Parse_Law_From_Args( char **argv, int index, int argc, Law_Desc The_Law )
{
  int		i;
  int		j;
  int		Nb_Params;
  boolean 	Error;
  boolean	Short;
  boolean 	Ok;
  double	Foo;
  
  i = index;
  Error = FALSE;
  Short = FALSE;

  The_Law->Name = argv[i][0];
  i++;

  if ( i < argc ) {
    switch( The_Law->Name ) {
    case 'I':
      /* this law has no parameters */
      break;
    case 'g':
    case 'G':
      if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[0] ) ) ) {
	Error = TRUE;
      }
      else
	i++;
      break;
    case 'A':
    case 'R':
    case 'W':
    case 'p':
    case 'N':
    case 'L':
    case 'l': /* double lognormale. Loi ad-hoc pour tests Experimental */
      The_Law->Parameters = (double*) 
	realloc ( The_Law->Parameters, 2*sizeof(double) );
      if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[0] ) ) ) {
	Error = TRUE;
      } 
      else {
	i++;
	if ( i < argc ) {
	  if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[1] ) ) ) {
	    Error = TRUE;
	  }
	  else 
	    i++;
	} 
	else {
	  Short = TRUE;
	}
      }
      break;
    case 'O':
      The_Law->Parameters = (double*) 
	realloc ( The_Law->Parameters, 3*sizeof(double) );
      if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[0] ) ) ) {
	Error = TRUE;
      }
      else {
	i++;
	if ( i < argc ) {
	  if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[1] ) ) ) {
	    Error = TRUE;
	  } 
	  else {
	    i++;
	    if ( i < argc ) {
	      if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[2] ) ) ) {
		Error = TRUE;
	      } 
	    }
	  }
	}
	else 
	  Error = TRUE;
      }
      break;
    case 'H':
      /* the number of parameters is not known in advance... have to 	*/
      /* find it first.							*/
      Ok = TRUE;
      for ( j=i; Ok && ( j < argc ); j++ ) {
	Ok = ( 1 == sscanf( argv[j], "%lf", &Foo ) );
      }
      /* counting the number of floats read */
      if ( j < argc )
	j--;
      Nb_Params = (j-i);
      if ( Nb_Params %2 != 0 ) {
	Error = TRUE;
      } else {
	The_Law->Parameters = (double*) 
	  realloc ( The_Law->Parameters, (1+Nb_Params)*sizeof(double) );
	The_Law->Parameters[0] = (double)(Nb_Params);
	for( j=i; j<i+Nb_Params; j++ ) {
	  sscanf( argv[j], "%lf", &(The_Law->Parameters[1+j-i]) );
	}
	i = i+Nb_Params;
      }
      break;
    case 'D':
      /* skip the (a priori) "[" argument */
      i++;
      /* the number of parameters is not known in advance... have to 	*/
      /* find it first, looking for the next "]".			*/
      Ok = TRUE;
      for ( j=i; Ok && ( j < argc ); j++ ) {
	Ok = ( 0 != strcmp( argv[j], "]" ) );
      }
      /* counting the number of floats read */
      if ( j >= argc ) {
	Short = TRUE;
      }
      else {
	j--;
	Nb_Params = (j-i);
	The_Law->Parameters = (double*) 
	  realloc ( The_Law->Parameters, (1+2*Nb_Params)*sizeof(double) );
	The_Law->Parameters[0] = (double)(Nb_Params);
	for( j=0; j<Nb_Params; j++ ) {
	  sscanf( argv[i+j], "%lf", &(The_Law->Parameters[1+j]) );
	}
	/* skip "]" and "[" */
	i += Nb_Params+2;
	/* continue reading. Have to protect against short arg list */
	Ok = TRUE;
	for( j=0; Ok && (j<Nb_Params) && (i+j<argc); j++ ) {
	  Ok = ( 1 == sscanf( argv[i+j], "%lf", 
			      &(The_Law->Parameters[1+Nb_Params+j]) ) );
	}
	if ( !Ok ) {
	  Error = TRUE;
	  i = i+j-1;
	} 
	else if ( i+j >= argc ) {
	  Short = TRUE;
	}
	else {
	  /* skip final "]", hoping it is actually there */
	  i = i + Nb_Params + 1;
	}
      }
      break;
    case 'P':
      /* read first parameter */
      if ( 1 != sscanf( argv[i], "%lf", &(The_Law->Parameters[0] ) ) ) {
	Error = TRUE;
      }
      else {
	i++;
	/* recursively read two laws */
	The_Law->Law_Parameters = (Law_Desc*) calloc( 2, sizeof(Law_Desc) );
	The_Law->Law_Parameters[0] = New_Law();
	The_Law->Law_Parameters[1] = New_Law();
	i = Parse_Law_From_Args( argv, i, argc, The_Law->Law_Parameters[0] );
	i = Parse_Law_From_Args( argv, i, argc, The_Law->Law_Parameters[1] );
      }
    default:
      printf("%s: unknown law id (%c)\n", TOOLNAME, The_Law->Name );
    }
  }
  else {
    Short = TRUE;
  }

  if ( Short ) {
    printf("%s: error parsing law %c: not enough arguments.\n", 
	   TOOLNAME, The_Law->Name );
    i = index;
  }
  else if ( Error ) {
    printf("%s: error parsing law %c at argument %s.\n", 
	   TOOLNAME, The_Law->Name, argv[i] );
    i = index;
  }

  return i;
}

Law_List Append_Law( Law_List	The_List, Law_Desc New_Law )
{
  Law_List Res;

  Res = (Law_List) malloc( sizeof( struct LAW_LIST ) );

  if ( Res != NULL ) {
    Res->Val = New_Law;
    Res->Next = The_List;
  }
  else {
    E_MEM( "Law list" );
  }

  return Res;
}

void Free_Law( Law_Desc The_Law )
{
  if ( The_Law->Parameters != NULL )
    free( The_Law->Parameters );
  
  /* number of Law parameters depends on the type */
  if ( The_Law->Name == 'P' ) {
    Free_Law( The_Law->Law_Parameters[0] );
    Free_Law( The_Law->Law_Parameters[1] );
    free( The_Law->Law_Parameters );
  }
    
  free( The_Law );
}
