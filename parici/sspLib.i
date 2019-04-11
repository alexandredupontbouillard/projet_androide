%module sspLib
%{

#include "header/marmoteSet.h"
#include "header/marmoteConstants.h"
#include "header/solutionMDP.h"
#include "header/Set/marmoteInterval.h"
#include "header/genericMDP.h"
#include "header/averageMDP.h"
#include "header/totalRewardMDP.h"
#include "header/transitionStructure.h"
#include "header/sparseMatrix.h"
#include "header/discountedMDP.h"
#include "header/feedbackSolutionMDP.h"
#include "header/nonStationarySolutionMDP.h"
#include "header/alglin.h"
%}

%include "header/marmoteSet.h"
%include "header/marmoteInterval.h"
%include "header/transitionStructure.h"
%include "header/sparseMatrix.h"
%include "header/marmoteConstants.h"
%include "header/solutionMDP.h"
%include "header/genericMDP.h"
%include "header/averageMDP.h"
%include "header/totalRewardMDP.h"
%include "std_string.i"
%include "header/alglin.h"
%include "header/nonStationarySolutionMDP.h"
%include "header/feedbackSolutionMDP.h"
%include "header/discountedMDP.h"
