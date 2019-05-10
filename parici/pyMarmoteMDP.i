%module pyMarmoteMDP 
%{
#include "header/marmoteSet.h"
#include "header/marmoteInterval.h"
#include "header/sparseMatrix.h"
#include "header/solutionMDP.h"
#include "header/feedbackSolutionMDP.h"
#include "header/genericMDP.h"
#include "header/totalRewardMDP.h"

%}

%include "std_string.i"
%include "std_vector.i"
%include "header/solutionMDP_SWIG.h"
%include "header/feedbackSolutionMDP_SWIG.h"
%include "header/marmoteSet_SWIG.h"
%include "header/marmoteInterval_SWIG.h"
%include "header/sparseMatrix_SWIG.h"
%include "header/totalRewardMDP_SWIG.h"

