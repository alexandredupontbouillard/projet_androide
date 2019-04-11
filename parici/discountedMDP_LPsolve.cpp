

#include <values.h>
#include "header/discountedMDP_LPsolve.h"
#include "alglin.h"

using namespace std;


// constructeurs

discountedMDP_LP::discountedMDP_LP(string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, sparseMatrix* rews, double b)
: discountedMDP(r,states,actions,trans,rews,b)
{

}


discountedMDP_LP::discountedMDP_LP(string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, vector<sparseMatrix*> rews, double b)
: discountedMDP(r,states,actions,trans,rews,b)
{


}

// destructeur
discountedMDP_LP::~discountedMDP_LP(){
}

// other utilities function
void discountedMDP_LP::writeMDP()
{
    discountedMDP::writeMDP();
    printf("#############################################\n");
    printf("Solving Using Linear Programming");
}


solutionMDP* linearProgrammingSolve(){ return NULL;}
