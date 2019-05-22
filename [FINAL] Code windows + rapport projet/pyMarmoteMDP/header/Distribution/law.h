/************************************************************************/
/* Headers for functions for manipulating distribution laws.
   This C structure is doomed to disappear, once all its functionalities
   will be reproduced in the marmoteCore Distribution object.		*/
/************************************************************************/
#ifndef _LAW_H
#define _LAW_H

#include <math.h>
#include "real_list.h"

/**
 * @brief Description of a distribution
 * @author Alain Jean-Maire
 *
 */
typedef struct LAW_DESC {
    char	Name; /**< name-identifier of the distribution */
    double*	Parameters; /**< array of standard parameters for the distribution */
    struct LAW_DESC	**Law_Parameters; /**< array of sub-distributions, e.g. for mixture or other compositions */
} *Law_Desc;

/**
 * @brief Chained list structure for lists of distributions
 * @author Alain Jean-Marie
 *
 */
typedef struct LAW_LIST {
    Law_Desc	Val; /**< the distribution at the no */
    struct LAW_LIST *Next; /**< pointer to the next in the list */
} *Law_List;


/* Definition of +oo for durations and rates */
#define INFINITE_DURATION       HUGE_VAL
#define INFINITE_RATE	        HUGE_VAL

/* Definitions for printing styles */
#define NO_PRINT_MODE		0
#define	DEFAULT_PRINT_MODE	1
#define	NORMAL_PRINT_MODE	1
#define QNAP_PRINT_MODE		2
#define PROSIT_PRINT_MODE	3
#define MMPP_PRINT_MODE		4
#define PNED_PRINT_MODE		5

/**
 * @brief
 *
 * @param Law_Desc
 * @return double
 */
double		Mean( Law_Desc );
/**
 * @brief
 *
 * @param Law_Desc
 * @return double
 */
double		Rate( Law_Desc );
/**
 * @brief
 *
 * @param Law_Desc
 * @param Order
 * @return double
 */
double		Moment( Law_Desc, int Order );
/**
 * @brief
 *
 * @param Law_Desc
 * @param s
 * @return double
 */
double		Laplace( Law_Desc, double s );
/**
 * @brief
 *
 * @param Law_Desc
 * @param s
 * @return double
 */
double		DLaplace( Law_Desc, double s );
/**
 * @brief
 *
 * @param Law_Desc
 * @param x
 * @return double
 */
double		CDF( Law_Desc, double x );
/**
 * @brief
 *
 * @param Law_Desc
 * @param Order
 * @return boolean
 */
boolean		hasMoment( Law_Desc, int Order );

/**
 * @brief
 *
 * @return Law_Desc
 */
Law_Desc New_Law( void );
/**
 * @brief
 *
 * @param Loi
 * @return Law_Desc
 */
Law_Desc Copy_Law( Law_Desc Loi );
/**
 * @brief
 *
 * @param Loi
 * @param Factor
 * @return Law_Desc
 */
Law_Desc Rescale_Law( Law_Desc Loi, double Factor );
/**
 * @brief
 *
 * @param Nom
 * @param Params
 * @param Params2
 * @param The_SubLaws
 * @return Law_Desc
 */
Law_Desc Parse_Law( char Nom, Liste_Reel Params,
		    Liste_Reel Params2, Law_List The_SubLaws );
/**
 * @brief
 *
 * @param
 * @param int
 * @param int
 * @param Law_Desc
 * @return int
 */
int Parse_Law_From_Args ( char**, int, int, Law_Desc );
/**
 * @brief
 *
 * @param Law_List
 * @param New_Law
 * @return Law_List
 */
Law_List Append_Law( Law_List, Law_Desc New_Law );
/**
 * @brief
 *
 * @param The_Law
 */
void Free_Law( Law_Desc The_Law );

/**
 * @brief
 *
 * @param Stream
 * @param Loi
 * @param Mode
 */
void	Write_Law( FILE *Stream, Law_Desc Loi, int Mode );
/* shortcut for printing, mostly out of backwards compatibility */
#define fprint_law(x)	Write_Law( stdout, x, NORMAL_PRINT_MODE )

#endif
